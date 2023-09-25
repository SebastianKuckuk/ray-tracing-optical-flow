#pragma once


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
