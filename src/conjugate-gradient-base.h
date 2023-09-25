#pragma once


//constexpr double gridWidth = 1. / 512.;
//constexpr double gridWidthSq = gridWidth * gridWidth;
//constexpr double gridWidthSqInv = 1. / gridWidthSq;
constexpr double dt = 1.e-0;
constexpr double regularization = 1.e4;//1.e4; // alpha**2 in the original formulation -> name clash with CG alpha


inline void initGradientsAndRHS(size_t nx, size_t ny, const double *__restrict__ img0, const double *__restrict__ img1,
                                double *__restrict__ ixix, double *__restrict__ ixiy, double *__restrict__ iyiy,
                                double *__restrict__ uRHS, double *__restrict__ vRHS) {

    auto gridWidth = 1. / (nx - 2.);

#pragma omp parallel for schedule (static) //collapse(2)
    for (size_t j = 1; j < ny - 1; ++j) {
        for (size_t i = 1; i < nx - 1; ++i) {
            auto ix = (img0[(j + 0) * nx + (i + 1)] - img0[j * nx + i]) / gridWidth;// + (img1[(j + 0) * nx + (i + 1)] - img1[j * nx + i]) / gridWidth;
            auto iy = (img0[(j + 1) * nx + (i + 0)] - img0[j * nx + i]) / gridWidth;// + (img1[(j + 1) * nx + (i + 0)] - img1[j * nx + i]) / gridWidth;
            auto it = (img1[j * nx + i] - img0[j * nx + i]) / dt;
            ixix[j * nx + i] = 1 * ix * ix;
            ixiy[j * nx + i] = 1 * ix * iy;
            iyiy[j * nx + i] = 1 * iy * iy;
            uRHS[j * nx + i] = -ix * it;
            vRHS[j * nx + i] = -iy * it;
        }
    }
}


inline void applyBC(size_t nx, size_t ny, double *__restrict__ u, double *__restrict__ v) {
    auto scale = -1.;

    // east
    for (size_t j = 1; j < ny - 1; ++j) {
        auto i = nx - 1;
        u[j * nx + i] = scale * u[(j + 0) * nx + (i - 1)];
        v[j * nx + i] = scale * v[(j + 0) * nx + (i - 1)];
    }

    // west
    for (size_t j = 1; j < ny - 1; ++j) {
        auto i = 0;
        u[j * nx + i] = scale * u[(j + 0) * nx + (i + 1)];
        v[j * nx + i] = scale * v[(j + 0) * nx + (i + 1)];
    }

    // north
    for (size_t i = 1; i < nx - 1; ++i) {
        auto j = ny - 1;
        u[j * nx + i] = scale * u[(j - 1) * nx + (i + 0)];
        v[j * nx + i] = scale * v[(j - 1) * nx + (i + 0)];
    }

    // south
    for (size_t i = 1; i < nx - 1; ++i) {
        auto j = 0;
        u[j * nx + i] = scale * u[(j + 1) * nx + (i + 0)];
        v[j * nx + i] = scale * v[(j + 1) * nx + (i + 0)];
    }
}


//inline void createImage(size_t nx, size_t ny,
//                   const double *const __restrict__ u, const double *const __restrict__ v,
//                   double *const __restrict__ uNew, double *const __restrict__ vNew,
//                   const double *const __restrict__ ixix, const double *const __restrict__ ixiy, const double *const __restrict__ iyiy,
//                   const double *const __restrict__ uRHS, const double *const __restrict__ vRHS) {
//    for (size_t j = 1; j < ny - 1; ++j) {
//        for (size_t i = 1; i < nx - 1; ++i) {
//            uNew[j * nx + i] = (uRHS[j * nx + i] - ixiy[j * nx + i] * v[j * nx + i]
//                                + regularization * ((u[j * nx + i - 1] + u[j * nx + i + 1] + u[(j - 1) * nx + i] + u[(j + 1) * nx + i]) / gridWidthSq))
//                               / (ixix[j * nx + i] + regularization * 4 / gridWidthSq);
//
//            vNew[j * nx + i] = (vRHS[j * nx + i] - ixiy[j * nx + i] * u[j * nx + i]
//                                + regularization * ((v[j * nx + i - 1] + v[j * nx + i + 1] + v[(j - 1) * nx + i] + v[(j + 1) * nx + i]) / gridWidthSq))
//                               / (iyiy[j * nx + i] + regularization * 4 / gridWidthSq);
//        }
//    }
//}


/*inline*/ size_t cg(size_t nx, size_t ny, size_t maxIt,
                     double *const __restrict__ u, double *const __restrict__ v,
                     double *const __restrict__ uRes, double *const __restrict__ vRes,
                     double *const __restrict__ uP, double *const __restrict__ vP,
                     double *const __restrict__ uAP, double *const __restrict__ vAP,
                     const double *const __restrict__ ixix, const double *const __restrict__ ixiy, const double *const __restrict__ iyiy,
                     const double *const __restrict__ uRHS, const double *const __restrict__ vRHS) {

    auto gridWidthSqInv = ((nx - 2.) * (nx - 2.));

    // initialization
    double initResSq = 0.;

#pragma omp parallel for schedule (static) reduction ( + : initResSq ) //collapse(2)
    for (size_t j = 1; j < ny - 1; ++j) {
        for (size_t i = 1; i < nx - 1; ++i) {
            // compute and store residual
            uRes[j * nx + i] = uRHS[j * nx + i] - (
                    ixix[j * nx + i] * u[j * nx + i]
                    + ixiy[j * nx + i] * v[j * nx + i]
                    + regularization * gridWidthSqInv *
                      (4 * u[j * nx + i] - (u[j * nx + i - 1] + u[j * nx + i + 1] + u[(j - 1) * nx + i] + u[(j + 1) * nx + i])));

            vRes[j * nx + i] = vRHS[j * nx + i] - (
                    ixiy[j * nx + i] * u[j * nx + i]
                    + iyiy[j * nx + i] * v[j * nx + i]
                    + regularization * gridWidthSqInv *
                      (4 * v[j * nx + i] - (v[j * nx + i - 1] + v[j * nx + i + 1] + v[(j - 1) * nx + i] + v[(j + 1) * nx + i])));

            // accumulate residual
            initResSq += uRes[j * nx + i] * uRes[j * nx + i] + vRes[j * nx + i] * vRes[j * nx + i];

            // init auxiliary fields
            uP[j * nx + i] = uRes[j * nx + i];
            vP[j * nx + i] = vRes[j * nx + i];
        }
    }
    applyBC(nx, ny, uRes, vRes);
    applyBC(nx, ny, uP, vP);

    // main loop
    double curResSq = initResSq;

    std::cout << 0 << '\t' << std::sqrt(curResSq) << std::endl;

    for (size_t it = 0; it < maxIt; ++it) {
        // compute A * p and alpha
        double alphaNominator = curResSq;
        double alphaDenominator = 0.;

#pragma omp parallel for schedule (static) reduction ( + : alphaDenominator ) //collapse(2)
        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t i = 1; i < nx - 1; ++i) {
                uAP[j * nx + i] = (
                        ixix[j * nx + i] * uP[j * nx + i]
                        + ixiy[j * nx + i] * vP[j * nx + i]
                        + regularization * gridWidthSqInv *
                          (4 * uP[j * nx + i] - (uP[j * nx + i - 1] + uP[j * nx + i + 1] + uP[(j - 1) * nx + i] + uP[(j + 1) * nx + i])));

                vAP[j * nx + i] = (
                        ixiy[j * nx + i] * uP[j * nx + i]
                        + iyiy[j * nx + i] * vP[j * nx + i]
                        + regularization * gridWidthSqInv *
                          (4 * vP[j * nx + i] - (vP[j * nx + i - 1] + vP[j * nx + i + 1] + vP[(j - 1) * nx + i] + vP[(j + 1) * nx + i])));

                alphaDenominator += uP[j * nx + i] * uAP[j * nx + i] + vP[j * nx + i] * vAP[j * nx + i];
            }
        }
        double alpha = alphaNominator / alphaDenominator;

        // update solution and residual
        double nextResSq = 0.;
#pragma omp parallel for schedule (static) reduction ( + : nextResSq ) //collapse(2)
        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t i = 1; i < nx - 1; ++i) {
                u[j * nx + i] += alpha * uP[j * nx + i];
                v[j * nx + i] += alpha * vP[j * nx + i];

                uRes[j * nx + i] -= alpha * uAP[j * nx + i];
                vRes[j * nx + i] -= alpha * vAP[j * nx + i];

                nextResSq += uRes[j * nx + i] * uRes[j * nx + i] + vRes[j * nx + i] * vRes[j * nx + i];
            }
        }
//        applyBC(nx, ny, u, v);
//        applyBC(nx, ny, uRes, vRes);

        // check exit criterion
        if (sqrt(nextResSq) <= 1e-8 * sqrt(initResSq)) {
            std::cout << it << '\t' << std::sqrt(nextResSq) << " <- " << std::sqrt(initResSq) << std::endl;
            return it;
        }

        // compute beta
        double beta = nextResSq / curResSq;
        curResSq = nextResSq;

        if (0 == it % 100)
            std::cout << it << '\t' << std::sqrt(curResSq) << std::endl;

        // update p
#pragma omp parallel for schedule (static) //collapse(2)
        for (size_t j = 1; j < ny - 1; ++j) {
            for (size_t i = 1; i < nx - 1; ++i) {
                uP[j * nx + i] = uRes[j * nx + i] + beta * uP[j * nx + i];
                vP[j * nx + i] = vRes[j * nx + i] + beta * vP[j * nx + i];
            }
        }
        applyBC(nx, ny, uP, vP);
    }

    std::cout << maxIt << '\t' << std::sqrt(curResSq) << std::endl;

    return maxIt;
}
