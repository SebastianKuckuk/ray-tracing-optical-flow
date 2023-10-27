#include <chrono>
#include <cstring>
#include <sstream>

#include <mpi.h>

#include "ray-tracing.h"
#include "optical-flow.h"
#include "image.h"
#include "timer.h"


void parseCLA(int argc, char *const *argv, size_t &numLevels, size_t &supersampling, size_t &mgIterations,
              double &tStart, double &dt, size_t &numTimeSteps, size_t &numReps, bool &printImages) {
    // default values
    numLevels = 10;
    supersampling = 2;
    mgIterations = 10;
    dt = 1e-2;
    numTimeSteps = 4;
    numReps = 1;
    printImages = false;

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
    if (argc > i) numTimeSteps = atoi(argv[i]);
    ++i;
    if (argc > i) numReps = atoi(argv[i]);
    ++i;
    if (argc > i) printImages = atoi(argv[i]);
    ++i;


    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    tStart = mpiRank * dt * numTimeSteps;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_INIT;
#endif

    size_t numLevels, supersampling, mgIterations;
    double tStart, dt;
    size_t numTimeSteps, numReps;
    bool printImages;
    parseCLA(argc, argv, numLevels, supersampling, mgIterations, tStart, dt, numTimeSteps, numReps, printImages);

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
    numSpheres = 48;
    spheres = static_cast<Sphere *>(malloc(numSpheres * sizeof(Sphere)));
    initSpheres();
    MPI_Bcast(spheres, numSpheres * (sizeof(Sphere) / sizeof(double)), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // prepare timers
    Timer timerRT("ray-tracing"), timerMap("map"), timerOF("optical-flow"), timerMPISend("mpi-send"), timerMPIRecv("mpi-recv");

    for (size_t rep = 0; rep < numReps; ++rep) {
        for (size_t tIt = 0; tIt < numTimeSteps + 1; ++tIt) {
            auto t = tStart + tIt * dt + rep * dt * numTimeSteps * numRanks;

            updateSpherePositions(t);

            // measurement
            if (numTimeSteps == tIt && !(numReps - 1 == rep && numRanks - 1 == mpiRank)) {
                timerMPIRecv.start();
                {
                    MPI_Recv(img, nxRT * nyRT, MPI_UNSIGNED_CHAR, (mpiRank + numRanks + 1) % numRanks, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                timerMPIRecv.stop();
            } else {
                timerRT.start();
                {
                    createImage(nxRT, nyRT, t, img);
                }
                timerRT.stop();
            }

            if (0 == tIt && !(0 == mpiRank && 0 == rep)) {
                timerMPISend.start();
                {
                    memcpy(mpiBuffer, img, nxRT * nyRT * sizeof(unsigned char) * sizeof(Color) / sizeof(double));
                    MPI_Isend(mpiBuffer, nxRT * nyRT, MPI_UNSIGNED_CHAR, (mpiRank + numRanks - 1) % numRanks, 0, MPI_COMM_WORLD, &mpiReq);
                }
                timerMPISend.stop();
            }

            timerMap.start();
            {
                mapImages(nxOF, nyOF, nxRT, img, img1, supersampling);
                std::swap(img0, img1);
            }
            timerMap.stop();

            if (tIt > 0) {
                timerOF.start();
                {
                    initGradientsAndRHS(nxOF, nyOF, dt, img0, img1, mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1], mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
                    for (auto level = numLevels - 1; level >= 1; --level)
                        coarsenOperator((1u << (level - 1)) + 2, (1u << (level - 1)) + 2, (1u << level) + 2, (1u << level) + 2,
                                        mgIxIx[level], mgIxIy[level], mgIyIy[level], mgIxIx[level - 1], mgIxIy[level - 1], mgIyIy[level - 1]);

//                auto initRes = resNorm(nxOF, nyOF,
//                                       mgSolU[numLevels - 1], mgSolV[numLevels - 1],
//                                       mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1],
//                                       mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);

                    for (auto mgIt = 0; mgIt < mgIterations; ++mgIt)
                        multigrid(numLevels - 1, 2. / 3., mgRhsU, mgRhsV, mgSolU, mgSolV, mgSolNewU, mgSolNewV, mgResU, mgResV, mgIxIx, mgIxIy, mgIyIy);

//                auto res = resNorm(nxOF, nyOF,
//                                   mgSolU[numLevels - 1], mgSolV[numLevels - 1],
//                                   mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1],
//                                   mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
//                std::cout << "Residual before / after " << mgIterations << " : " << initRes << " / " << res << std::endl;
                }
                timerOF.stop();
            }

            // print images
            if (printImages) {
                if (tIt < numTimeSteps) {
                    std::stringstream filename;
                    filename << "../images/ray-tracing-" << std::setw(5) << std::setfill('0') << std::right << t / dt << ".raw";
                    printImage(nxOF, nyOF, img0, filename.str());
                }
                if (tIt > 0) {
                    std::stringstream filename;
                    filename << "../images/optical-flow-" << std::setw(5) << std::setfill('0') << std::right << t / dt << ".raw";
                    printImage(nxOF, nyOF, mgSolU[numLevels - 1], mgSolV[numLevels - 1], filename.str());
                }
            }
        }

        if (mpiRank > 0 || rep > 0) {
            MPI_Wait(&mpiReq, MPI_STATUS_IGNORE);
            mpiReq = MPI_REQUEST_NULL;
        }
    }

    for (auto timer: {timerRT, timerMap, timerOF, timerMPISend, timerMPIRecv}) {
        for (auto r = 0; r < numRanks; ++r) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (mpiRank == r)
                timer.print(mpiRank);
            MPI_Barrier(MPI_COMM_WORLD); // TODO: remove
        }
    }

    // de-allocate MPI
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

#ifdef LIKWID_PERFMON
    LIKWID_MARKER_CLOSE;
#endif

    return 0;
}
