#include <chrono>
#include <fstream>
#include <cstring>
#include <sstream>

#include <mpi.h>

#include "ray-tracing.h"
#include "optical-flow.h"


void parseCLA(int argc, char *const *argv, size_t &numLevels, size_t &supersampling, size_t &mgIterations,
              double &tStart, double &dt, size_t &numTimeSteps) {
    // default values
    numLevels = 10;
    supersampling = 2;
    mgIterations = 10;
    dt = 1e-2;
    numTimeSteps = 4;

    // override with command line arguments
    int i = 1;
    if (argc > i) numLevels = atoi(argv[i]);
    ++i;
    if (argc > i) supersampling = atoi(argv[i]);
    ++i;
    if (argc > i) mgIterations = atoi(argv[i]);
    ++i;
    if (argc > i) dt = atof(argv[i]);
    ++i;
    if (argc > i) numTimeSteps = atof(argv[i]);
    ++i;

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    tStart = mpiRank * dt * numTimeSteps;
}


void mapImages(size_t nxOF, size_t nyOF, size_t nxRT, const unsigned char *img, double *img0, size_t supersampling) {
#pragma omp parallel for schedule (static) collapse(2)
    for (size_t j = 0; j < nyOF; ++j)
        for (size_t i = 0; i < nxOF; ++i) {
            img0[j * nxOF + i] = 0;

            for (size_t jOff = 0; jOff < supersampling; ++jOff) {
                for (size_t iOff = 0; iOff < supersampling; ++iOff) {
#ifdef USE_COLOR
                    for (auto dim = 0; dim < 3; ++dim)
                        img0[j * nxOF + i] += img[((supersampling * j + jOff) * nxRT + supersampling * i + iOff) * 3 + dim];
                    img0[j * nxOF + i] /= 3;
#else
                    img0[j * nxOF + i] += img[(supersampling * j + jOff) * nxRT + supersampling * i + iOff];
#endif
                }
            }

            img0[j * nxOF + i] /= supersampling * supersampling;
            img0[j * nxOF + i] /= 255.;
        }
}


void printImage(size_t nx, size_t ny, const double *const img, const std::string &filename) {
    std::ofstream outStream(filename, std::iostream::binary);
    for (size_t j = 1; j < ny - 1; ++j) {
        for (size_t i = 1; i < nx - 1; ++i) {
            for (auto dim = 0; dim < 3; ++dim) {
                auto c = (float) img[j * nx + i];
                outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            }
        }
    }

    outStream.close();
}

void printImage(size_t nx, size_t ny, const double *const imgU, const double *const imgV, const std::string &filename) {
    std::ofstream outStream(filename, std::iostream::binary);
    for (size_t j = 1; j < ny - 1; ++j) {
        for (size_t i = 1; i < nx - 1; ++i) {
            float c;
            c = (float) (1. * imgU[j * nx + i] + 0.);
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            c = (float) (1. * imgV[j * nx + i] + 0.);
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            c = (float) (0.);
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
        }
    }

    outStream.close();
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    size_t numLevels, supersampling, mgIterations;
    double tStart, dt;
    size_t numTimeSteps;
    parseCLA(argc, argv, numLevels, supersampling, mgIterations, tStart, dt, numTimeSteps);

    size_t nxOF = (1u << (numLevels - 1)) + 2;
    size_t nyOF = (1u << (numLevels - 1)) + 2;
    auto nxRT = supersampling * nxOF;
    auto nyRT = supersampling * nyOF;

    // allocate for ray tracer
    auto img = new unsigned char[nxRT * nyRT * sizeof(Color) / sizeof(double)];

    // allocate for optical flow solver
    auto img0 = new double[nxOF * nyOF];
    auto img1 = new double[nxOF * nyOF];

    auto mgSolU = new double *[numLevels];
    auto mgSolV = new double *[numLevels];
    auto mgSolNewU = new double *[numLevels];
    auto mgSolNewV = new double *[numLevels];
    auto mgResU = new double *[numLevels];
    auto mgResV = new double *[numLevels];
    auto mgRhsU = new double *[numLevels];
    auto mgRhsV = new double *[numLevels];
    auto mgIxIx = new double *[numLevels];
    auto mgIxIy = new double *[numLevels];
    auto mgIyIy = new double *[numLevels];

    for (auto level = 0; level < numLevels; ++level) {
        auto nx = (1u << level) + 2;
        auto ny = (1u << level) + 2;

        for (auto ptr: {mgSolU, mgSolV, mgSolNewU, mgSolNewV, mgResU, mgResV, mgRhsU, mgRhsV, mgIxIx, mgIxIy, mgIyIy}) {
            ptr[level] = new double[nx * ny];
            memset(ptr[level], 0, nx * ny);
        }
    }

    // init MPI
    int mpiRank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    MPI_Request mpiReq;
    auto mpiBuffer = new unsigned char[nxRT * nyRT * sizeof(Color) / sizeof(double)];

    // allocate spheres
    numSpheres = 1 + 1 + 20; // one large bottom sphere, one center sphere and multiple orbiting spheres
    spheres = static_cast<Sphere *>(malloc(numSpheres * sizeof(Sphere)));

    for (size_t tIt = 0; tIt < numTimeSteps + 1; ++tIt) {
        auto t = tStart + tIt * dt;

        initSpheres(t);

        // measurement
        auto startRT = std::chrono::steady_clock::now();

        if (numTimeSteps == tIt && mpiRank < numRanks - 1)
            MPI_Recv(img, nxRT * nyRT, MPI_UINT8_T, mpiRank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else
            createImage(nxRT, nyRT, t, img);

        if (0 == tIt && mpiRank > 0) {
            memcpy(mpiBuffer, img, nxRT * nyRT * sizeof(unsigned char) * sizeof(Color) / sizeof(double));
            MPI_Isend(mpiBuffer, nxRT * nyRT, MPI_UINT8_T, mpiRank - 1, 0, MPI_COMM_WORLD, &mpiReq);
        }

        auto endRT = std::chrono::steady_clock::now();

        auto startMap = std::chrono::steady_clock::now();

        std::swap(img0, img1);

        mapImages(nxOF, nyOF, nxRT, img, img0, supersampling);

        auto endMap = std::chrono::steady_clock::now();

        auto startOF = std::chrono::steady_clock::now();

        if (tIt > 0) {
            initGradientsAndRHS(nxOF, nyOF, dt, img0, img1, mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1], mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
            for (auto level = numLevels - 1; level >= 1; --level)
                coarsenOperator((1u << (level - 1)) + 2, (1u << (level - 1)) + 2, (1u << level) + 2, (1u << level) + 2,
                                mgIxIx[level], mgIxIy[level], mgIyIy[level], mgIxIx[level - 1], mgIxIy[level - 1], mgIyIy[level - 1]);

            auto initRes = resNorm(nxOF, nyOF,
                                   mgSolU[numLevels - 1], mgSolV[numLevels - 1],
                                   mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1],
                                   mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);

            for (auto mgIt = 0; mgIt < mgIterations; ++mgIt)
                multigrid(numLevels - 1, 2. / 3., mgRhsU, mgRhsV, mgSolU, mgSolV, mgSolNewU, mgSolNewV, mgResU, mgResV, mgIxIx, mgIxIy, mgIyIy);

            auto res = resNorm(nxOF, nyOF,
                               mgSolU[numLevels - 1], mgSolV[numLevels - 1],
                               mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1],
                               mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
            std::cout << "Residual before / after " << mgIterations << " : " << initRes << " / " << res << std::endl;
        }

        auto endOF = std::chrono::steady_clock::now();

        std::chrono::duration<double> elapsedSecondsRT = endRT - startRT;
        std::chrono::duration<double> elapsedSecondsMap = endMap - startMap;
        std::chrono::duration<double> elapsedSecondsOF = endOF - startOF;
        std::cout << "Time for ray-tracing / mapping / optical flow "
                  << 1e3 * elapsedSecondsRT.count() << " / "
                  << 1e3 * elapsedSecondsMap.count() << " / "
                  << 1e3 * elapsedSecondsOF.count() << std::endl;

        // print images
        if (tIt < numTimeSteps) {
            std::stringstream filename;
            filename << "../images/ray-tracing-" << t << ".raw";
            printImage(nxOF, nyOF, img0, filename.str());
        }
        if (tIt > 0) {
            std::stringstream filename;
            filename << "../images/optical-flow-" << t - 0.5 * dt << ".raw";
            printImage(nxOF, nyOF, mgSolU[numLevels - 1], mgSolV[numLevels - 1], filename.str());
        }
    }

    // de-allocate MPI
    if (mpiRank > 0)
        MPI_Wait(&mpiReq, MPI_STATUS_IGNORE);
    delete[] mpiBuffer;

    // de-allocate optical flow solver
    delete[] img0;
    delete[] img1;

    for (auto ptr: {mgSolU, mgSolV, mgSolNewU, mgSolNewV,
                    mgResU, mgResV, mgRhsU, mgRhsV, mgIxIx,
                    mgIxIy, mgIyIy}) {
        for (auto level = 0; level < numLevels; ++level)
            delete[] ptr[level];
        delete[] ptr;
    }

    // de-allocate ray tracer
    delete[] img;
    free(spheres);

    MPI_Finalize();

    return 0;
}
