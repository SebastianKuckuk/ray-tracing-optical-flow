#pragma once

#include <fstream>
#include <cstring>


void mapImageToDouble(size_t nxOF, size_t nyOF, size_t nxRT, const Color *imgSrc, double *imgDest, size_t supersampling) {
#pragma omp parallel for schedule (static) collapse(2)
    for (size_t j = 0; j < nyOF; ++j)
        for (size_t i = 0; i < nxOF; ++i) {
            imgDest[j * nxOF + i] = 0;

            for (size_t jOff = 0; jOff < supersampling; ++jOff) {
                for (size_t iOff = 0; iOff < supersampling; ++iOff) {
#ifdef USE_COLOR
                    imgDest[j * nxOF + i] += 1. / 3. * (imgSrc[(supersampling * j + jOff) * nxRT + supersampling * i + iOff].x
                                                        + imgSrc[(supersampling * j + jOff) * nxRT + supersampling * i + iOff].y
                                                        + imgSrc[(supersampling * j + jOff) * nxRT + supersampling * i + iOff].z);
#else
                    imgDest[j * nxOF + i] += imgSrc[(supersampling * j + jOff) * nxRT + supersampling * i + iOff];
#endif
                }
            }

            imgDest[j * nxOF + i] /= supersampling * supersampling;
        }
}


void mapImageToColor(size_t nxOF, size_t nyOF, size_t nxRT, const Color *imgSrc, Color *imgDest, size_t supersampling) {
#pragma omp parallel for schedule (static) collapse(2)
    for (size_t j = 0; j < nyOF; ++j) {
        for (size_t i = 0; i < nxOF; ++i) {
            imgDest[j * nxOF + i] = Color{0};

            for (size_t jOff = 0; jOff < supersampling; ++jOff)
                for (size_t iOff = 0; iOff < supersampling; ++iOff)
                    imgDest[j * nxOF + i] += imgSrc[(supersampling * j + jOff) * nxRT + supersampling * i + iOff];

            imgDest[j * nxOF + i] = imgDest[j * nxOF + i] / (supersampling * supersampling);
        }
    }
}


void printImage(size_t nx, size_t ny, const Color *const img, const std::string &filename) {
    std::ofstream outStream(filename, std::iostream::binary);
    for (size_t j = 1; j < ny - 1; ++j) {
        for (size_t i = 1; i < nx - 1; ++i) {
#ifdef USE_COLOR
            auto c = (float) img[j * nx + i].x;
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            c = (float) img[j * nx + i].y;
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            c = (float) img[j * nx + i].z;
            outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
#else
            for (auto dim = 0; dim < 3; ++dim) {
                auto c = (float) img[j * nx + i];
                outStream.write(reinterpret_cast<const char *>(&c), sizeof(float));
            }
#endif
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
