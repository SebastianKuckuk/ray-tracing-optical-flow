#include <chrono>
#include <fstream>
#include <cstring>
#include <sstream>

#include "ray-tracing.h"
#include "optical-flow.h"


void parseCLA(int argc, char *const *argv, size_t &numLevels, size_t &supersampling, size_t &mgIterations, double &dt, double &maxTime) {
    // default values
    numLevels = 10;
    supersampling = 2;
    mgIterations = 10;
    dt = 1e-2;
    maxTime = dt * 2;

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
    if (argc > i) maxTime = atof(argv[i]);
    ++i;
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
            c = (float) (2. * imgU[j * nx + i] + 0.5);
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            c = (float) (2. * imgV[j * nx + i] + 0.5);
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            c = (float) (0.);
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
        }
    }

    outStream.close();
}


int main(int argc, char *argv[]) {
    size_t numLevels, supersampling, mgIterations;
    double dt, maxTime;
    parseCLA(argc, argv, numLevels, supersampling, mgIterations, dt, maxTime);

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

    auto t = 0.0;

    // init spheres
    numSpheres = 1 + 1 + 20; // one large bottom sphere, one center sphere and multiple orbiting spheres
    spheres = static_cast<Sphere *>(malloc(numSpheres * sizeof(Sphere)));
    spheres[0] = {{0., 0., -130.}, 128 * 128., Color{0.8}, Color{0.2}, 16, Color{0.8}};
    spheres[1] = {{0., 0., 0.}, 0.5 * 0.5, Color{0.8}, Color{0.2}, 8, Color{0.2}};
    // remaining spheres will be initialized in the time loop since they are moving

    for (; t <= maxTime; t += dt) {
        for (auto i = 2; i < numSpheres; ++i) {
            auto r = 1. / 8.;

            auto alpha = 2. * M_PI * (-(t + (double) (i - 2) / (numSpheres - 2)));
            auto beta = 2. * M_PI * ((i - 2 >= (numSpheres - 2) / 2 ? -1 : 1) / 16.);
            auto gamma = 2. * M_PI * (0.);

            auto pos = Vec3{sqrt(spheres[1].rSq) + ((numSpheres - 2) / 4.) * r, 0, 0}.rotate(alpha, beta, gamma);

            spheres[i] = {pos, r * r, Color{0.2}, Color{0.8}, 32, Color{0.8}};
        }

        // measurement
        auto startRT = std::chrono::steady_clock::now();

        createImage(nxRT, nyRT, t, img);

        auto endRT = std::chrono::steady_clock::now();

        auto startMap = std::chrono::steady_clock::now();

        std::swap(img0, img1);

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

        auto endMap = std::chrono::steady_clock::now();

        auto startOF = std::chrono::steady_clock::now();

        if (t > 0) {
            initGradientsAndRHS(nxOF, nyOF, dt, img0, img1, mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1], mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
            for (auto level = numLevels - 1; level >= 1; --level)
                coarsenOperator((1u << (level - 1)) + 2, (1u << (level - 1)) + 2, (1u << level) + 2, (1u << level) + 2,
                                mgIxIx[level], mgIxIy[level], mgIyIy[level], mgIxIx[level - 1], mgIxIy[level - 1], mgIyIy[level - 1]);

            auto initRes = resNorm(nxOF, nyOF,
                                   mgSolU[numLevels - 1], mgSolV[numLevels - 1],
                                   mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1],
                                   mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
            std::cout << "Initial residual : " << initRes << std::endl;

            for (auto mgIt = 0; mgIt < mgIterations; ++mgIt)
                multigrid(numLevels - 1, 2. / 3., mgRhsU, mgRhsV, mgSolU, mgSolV, mgSolNewU, mgSolNewV, mgResU, mgResV, mgIxIx, mgIxIy, mgIyIy);

            auto res = resNorm(nxOF, nyOF,
                               mgSolU[numLevels - 1], mgSolV[numLevels - 1],
                               mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1],
                               mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
            std::cout << "Residual after " << mgIterations << " : " << res << std::endl;
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
        {
            std::stringstream filename;
            filename << "../images/ray-tracing-" << t << ".raw";
            printImage(nxOF, nyOF, img0, filename.str());
        }
        {
            std::stringstream filename;
            filename << "../images/optical-flow-" << t << ".raw";
            printImage(nxOF, nyOF, mgSolU[numLevels - 1], mgSolV[numLevels - 1], filename.str());
        }
    }

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

    return 0;
}
