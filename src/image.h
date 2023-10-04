#pragma once

#include <fstream>
#include <cstring>


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
    //for (size_t j = 1; j < ny - 1; ++j) {
    for (size_t j = ny - 2; j > 0; --j) {
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
    //for (size_t j = 1; j < ny - 1; ++j) {
    for (size_t j = ny - 2; j > 0; --j) {
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
