#include <chrono>
#include <fstream>
#include <cstring>

#define USE_CIMG
#ifdef USE_CIMG

#include "Cimg.h"

#endif

#include "ray-tracing-util.h"
#include "intersect.h"
#include "conjugate-gradient-base.h"
#include "multigrid-base.h"


Color trace(const Ray &ray, int numBounces) {
    // locate the closest intersection
    int iHit = 0;
    double rayDist = 1e6;
    if (!intersect(ray, rayDist, iHit))
        return Color{0};

    Sphere &sphere = spheres[iHit];

    auto hitPoint = ray.origin + rayDist * ray.direction;
    auto normal = hitPoint - sphere.pos;
    normal = normal.normalize();

    auto reflectionRay = Ray{hitPoint, ray.direction - 2 * dot(normal, ray.direction) * normal};

    // add lightning
    Color color(0);

    for (const auto &light: lights) {
        auto lightRay = Ray{hitPoint, (light.position - hitPoint).normalize()};

        // check if light is visible - first compare light direction with surface normal
        auto diffuseIntensity = dot(lightRay.direction, normal);
        if (diffuseIntensity < 0)
            continue;

        // then check if the light source is occluded
        auto lightRayDist = 1e6;
        if (intersect(lightRay, lightRayDist, iHit))
            continue;

        // add diffuse contribution
        color += diffuseIntensity * sphere.diffuseColor * light.color;

        // add specular contribution
        auto specularIntensity = dot(reflectionRay.direction, lightRay.direction);
        if (specularIntensity > 0)
            color += pow(specularIntensity, sphere.specularExp) * sphere.specularColor * light.color;
    }

    // add contribution from reflection ray
    if (numBounces < 8)
        color += sphere.reflectionColor * trace(reflectionRay, numBounces + 1);

    return color;
}


inline void createImage(size_t nx, size_t ny, double t, unsigned char *__restrict__ img) {
#pragma omp parallel for
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
//            auto x = (double) i / (double) (nx - 1);
//            auto y = 1. - (double) j / (double) (ny - 1);
//
//            auto direction = Vec3{-0.847569 - x * 1.30741 - y * 1.19745, -1.98535 + x * 2.11197 - y * 0.741279, -2.72303 + y * 2.04606};
//
//            Ray ray{{2.1, 1.3, 1.7}, direction.normalize()};
//            Color color = trace(ray, 0);

            auto x = (i + 0.5) / nx;
            auto y = (j + 0.5) / nx;

            auto fov = 90;
            auto scale = tan(fov * 0.5 * M_PI / 180);
            auto imageAspectRatio = nx / (double) ny;

            double alpha = t * 2. * M_PI, beta = 0., gamma = -0.25 * M_PI;
            double rot[3][3] = {
                    {
                            cos(alpha) * cos(beta),
                            cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma),
                            cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma)
                    },
                    {
                            sin(alpha) * cos(beta),
                            sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma),
                            sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma)
                    },
                    {
                            -sin(beta),
                            cos(beta) * sin(gamma),
                            cos(beta) * cos(gamma)
                    }
            };

            auto origin = Vec3{0., 0., 2.}; // TODO: camMat * vec
//            auto origin = Vec3{2.1, 1.3, 1.7};

            auto px = (2 * (i + 0.5) / nx - 1) * scale * imageAspectRatio;
            auto py = (1 - 2 * (j + 0.5) / ny) * scale;
            auto direction = Vec3{px, py, -1};// - origin;// TODO: camMat * vec

            origin = Vec3{
                    rot[0][0] * origin.x + rot[0][1] * origin.y + rot[0][2] * origin.z,
                    rot[1][0] * origin.x + rot[1][1] * origin.y + rot[1][2] * origin.z,
                    rot[2][0] * origin.x + rot[2][1] * origin.y + rot[2][2] * origin.z
            };

            direction = Vec3{
                    rot[0][0] * direction.x + rot[0][1] * direction.y + rot[0][2] * direction.z,
                    rot[1][0] * direction.x + rot[1][1] * direction.y + rot[1][2] * direction.z,
                    rot[2][0] * direction.x + rot[2][1] * direction.y + rot[2][2] * direction.z
            };

            Ray ray{origin, direction.normalize()};
            Color color = trace(ray, 0);

            // clamp color to [0, 1] and store result
#ifdef USE_COLOR
            color.x = std::min(1., std::max(0., color.x));
            color.y = std::min(1., std::max(0., color.y));
            color.z = std::min(1., std::max(0., color.z));
            img[(j * nx + i) * 3 + 0] = (unsigned char) (255 * color.x);
            img[(j * nx + i) * 3 + 1] = (unsigned char) (255 * color.y);
            img[(j * nx + i) * 3 + 2] = (unsigned char) (255 * color.z);
#else
            color = std::min(1., std::max(0., color));
            img[j * nx + i] = (unsigned char) (255 * color);
#endif
        }
    }
}


int main(int argc, char *argv[]) {
    size_t nxOF, nyOF, nItWarmUp, nIt;
    parseCLA_2d(argc, argv, nxOF, nyOF, nItWarmUp, nIt);

    auto nxRT = 2 * nxOF;
    auto nyRT = 2 * nyOF;

    // allocate for ray tracer
    auto img = new unsigned char[nxRT * nyRT * sizeof(Color) / sizeof(double)];

    // allocate for optical flow solver
    auto img0 = new double[nxOF * nyOF];
    auto img1 = new double[nxOF * nyOF];

    auto ixix = new double[nxOF * nyOF];
    auto ixiy = new double[nxOF * nyOF];
    auto iyiy = new double[nxOF * nyOF];

    auto u = new double[nxOF * nyOF];
    auto v = new double[nxOF * nyOF];
    auto uRHS = new double[nxOF * nyOF];
    auto vRHS = new double[nxOF * nyOF];
    auto uRes = new double[nxOF * nyOF];
    auto vRes = new double[nxOF * nyOF];

    auto uP = new double[nxOF * nyOF];
    auto vP = new double[nxOF * nyOF];
    auto uAP = new double[nxOF * nyOF];
    auto vAP = new double[nxOF * nyOF];

    // multigrid
    auto numLevels = 10;
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

        mgSolU[level] = new double[nx * ny];
        mgSolV[level] = new double[nx * ny];
        mgSolNewU[level] = new double[nx * ny];
        mgSolNewV[level] = new double[nx * ny];
        mgResU[level] = new double[nx * ny];
        mgResV[level] = new double[nx * ny];
        mgRhsU[level] = new double[nx * ny];
        mgRhsV[level] = new double[nx * ny];

        mgIxIx[level] = new double[nx * ny];
        mgIxIy[level] = new double[nx * ny];
        mgIyIy[level] = new double[nx * ny];

        for (auto ptr: {mgSolU[level], mgSolV[level], mgSolNewU[level], mgSolNewV[level],
                        mgResU[level], mgResV[level], mgRhsU[level], mgRhsV[level],
                        mgIxIx[level], mgIxIy[level], mgIyIy[level]})
            memset(ptr, 0, nx * ny);
    }

    // init
    memset(u, 0, nxOF * nxOF * sizeof(double));
    memset(v, 0, nxOF * nxOF * sizeof(double));

    auto t = 0.0;

    // init
    numSpheres = 2 + 20;
    spheres = static_cast<Sphere *>(malloc(numSpheres * sizeof(Sphere)));
    spheres[0] = {{0., 0., -130.}, 128 * 128., Color{0.8}, Color{0.2}, 16, Color{0.8}};
    spheres[1] = {{0., 0., 0.}, 0.5 * 0.5, Color{0.8}, Color{0.2}, 8, Color{0.2}};

#ifdef USE_CIMG
    for (; t <= 1000; t += 1.e-2) {
#else
        for (; t <= 1e-3; t += 1e-3) {
#endif
        for (auto i = 2; i < numSpheres; ++i) {
//        auto scale = 2. * M_PI / (numSpheres - 2);
//        auto pos = Vec3{cos(t + scale * (i - 2)),
//                        sin(t + scale * (i - 2)) * cos(0.1 * M_PI),
//                        sin(t + scale * (i - 2)) * sin(0.1 * M_PI)};

            auto alpha = (t + (double) (i - 2) / (numSpheres - 2)) * 2. * M_PI;
            auto beta = (i - 2 >= (numSpheres - 2) / 2 ? -0.125 : 0.125) * 2. * M_PI;
            auto gamma = 0. * 2. * M_PI;

            auto pos = Vec3{
                    cos(alpha) * cos(beta),
                    cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma),
                    cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma)
            };

            auto r = 1. / 6.;
            pos *= 0.5 + 1.5 * r;
            spheres[i] = {pos, r * r, Color{0.2}, Color{0.8}, 32, Color{0.8}};
        }

        // measurement
        auto start = std::chrono::steady_clock::now();

        auto startRT = std::chrono::steady_clock::now();

        createImage(nxRT, nyRT, t, img);
        std::swap(img0, img1);

        for (size_t j = 0; j < nyOF; ++j)
            for (size_t i = 0; i < nxOF; ++i)
#ifdef USE_COLOR
            {
                img0[j * nxOF + i] = 0;
                for (auto dim = 0; dim < 3; ++dim)
                    img0[j * nxOF + i] += img[((2 * j + 0) * nxRT + 2 * i + 0) * 3 + dim]
                                          + img[((2 * j + 0) * nxRT + 2 * i + 1) * 3 + dim]
                                          + img[((2 * j + 1) * nxRT + 2 * i + 0) * 3 + dim]
                                          + img[((2 * j + 1) * nxRT + 2 * i + 1) * 3 + dim];
                img0[j * nxOF + i] /= 3 * 4;
//                for (auto dim = 0; dim < 3; ++dim)
//                    img0[j * nxOF + i] += img[(j * nxRT + i) * 3 + dim];
//                img0[j * nxOF + i] /= 3;
            }
#else
        img0[j * nxOF + i] += img[(2 * j + 0) * nxRT + 2 * i + 0]
                              + img[(2 * j + 0) * nxRT + 2 * i + 1]
                              + img[(2 * j + 1) * nxRT + 2 * i + 0]
                              + img[(2 * j + 1) * nxRT + 2 * i + 1];
//        img0[j * nxOF + i] = img[j * nxOF + i];
#endif

        auto endRT = std::chrono::steady_clock::now();

        auto startOF = std::chrono::steady_clock::now();

        if (t > 0) {
            size_t nIt = 0;
            auto maxIt = 1024;
            applyBC(nxOF, nyOF, u, v);
            initGradientsAndRHS(nxOF, nyOF, img0, img1, ixix, ixiy, iyiy, uRHS, vRHS);
//            for (; nIt < maxIt;)

//            nIt += cg(nxOF, nyOF, maxIt, u, v, uRes, vRes, uP, vP, uAP, vAP, ixix, ixiy, iyiy, uRHS, vRHS);

            memcpy(mgSolU[numLevels - 1], u, nxOF * nyOF * sizeof(double));
            memcpy(mgSolV[numLevels - 1], v, nxOF * nyOF * sizeof(double));
            memcpy(mgRhsU[numLevels - 1], uRHS, nxOF * nyOF * sizeof(double));
            memcpy(mgRhsV[numLevels - 1], vRHS, nxOF * nyOF * sizeof(double));
            memcpy(mgIxIx[numLevels - 1], ixix, nxOF * nyOF * sizeof(double));
            memcpy(mgIxIy[numLevels - 1], ixiy, nxOF * nyOF * sizeof(double));
            memcpy(mgIyIy[numLevels - 1], iyiy, nxOF * nyOF * sizeof(double));
            for (auto mgIt = 0; mgIt < 10; ++mgIt)
                multigrid(numLevels - 1, 2. / 3., mgRhsU, mgRhsV, mgSolU, mgSolV, mgSolNewU, mgSolNewV, mgResU, mgResV, mgIxIx, mgIxIy, mgIyIy);
            memcpy(u, mgSolU[numLevels - 1], nxOF * nyOF * sizeof(double));
            memcpy(v, mgSolV[numLevels - 1], nxOF * nyOF * sizeof(double));

//            nIt += cg(nxOF, nyOF, maxIt, u, v, uRes, vRes, uP, vP, uAP, vAP, ixix, ixiy, iyiy, uRHS, vRHS);

            std::cout << "CG steps: " << nIt << std::endl;
        }

        auto endOF = std::chrono::steady_clock::now();

        std::chrono::duration<double> elapsedSecondsRT = endRT - startRT;
        std::chrono::duration<double> elapsedSecondsOF = endOF - startOF;
        std::cout << "Time for ray-tracing / optical flow " << 1e3 * elapsedSecondsRT.count() << " / " << 1e3 * elapsedSecondsOF.count() << std::endl;

        auto end = std::chrono::steady_clock::now();

#ifdef USE_CIMG
        cimg_library::CImg<double> imgShowRT(nxRT, nyRT, 1, 3, 1.);
        cimg_library::CImg<double> imgShowOF(nxOF, nyOF, 1, 3, 1.);
        for (size_t j = 0; j < nyRT; ++j)
            for (size_t i = 0; i < nxRT; ++i)
                for (auto dim = 0; dim < 3; ++dim)
#ifdef USE_COLOR
                        *imgShowRT.data(i, j, dim, 0) = img[(j * nxRT + i) * 3 + dim];
#else
        *imgShowRT.data(i, j, dim, 0) = img[j * nxRT + i];
#endif

        for (size_t j = 0; j < nyOF; ++j)
            for (size_t i = 0; i < nxOF; ++i)
                for (auto dim = 0; dim < 3; ++dim)
                    *imgShowOF.data(i, j, dim, 0) = std::sqrt(u[j * nxOF + i] * u[j * nxOF + i] + v[j * nxOF + i] * v[j * nxOF + i]);
        //                    *imgShowOF.data(i, j, dim, 0) = std::sqrt(uRes[j * nxOF + i] * uRes[j * nxOF + i] + vRes[j * nxOF + i] * vRes[j * nxOF + i]);

        imgShowRT.mirror('y');
        imgShowOF.mirror('y');
        imgShowRT.resize(512, 512);
        imgShowOF.resize(512, 512);
        static cimg_library::CImgDisplay displayRT(imgShowRT, "ray-tracing");
        displayRT = imgShowRT;
        static cimg_library::CImgDisplay displayOF(imgShowRT, "optical-flow");
        displayOF = imgShowOF;
        while (t >= 1000 && !displayRT.is_closed() && !displayOF.is_closed())
            displayRT.wait();
    }
#else
    }
//    printStats(end - start, nxOF * nyOF, nIt, stencil2dNumReads, stencil2dNumWrites);
#endif

    // check solution
#ifdef USE_COLOR
    std::ofstream outStream("./result.raw");
    for (size_t j = 0; j < nyOF; ++j)
        for (size_t i = 0; i < nxOF; ++i)
            outStream << img[(j * nxOF + i) * 3 + 0] << img[(j * nxOF + i) * 3 + 1] << img[(j * nxOF + i) * 3 + 2];
    outStream.close();
#else
    //    checkSolutionRayTracing(img, nxOF, nyOF);
    auto fd = fopen("./result.pnm", "w");
    fprintf(fd, "P5\n%zu %zu\n255\n", nxOF, nyOF);
    fwrite(img, sizeof(unsigned char), nxOF * nyOF, fd);
    fclose(fd);

    std::ofstream outStream("./result.raw");
    for (size_t j = 0; j < nyOF; ++j)
        for (size_t i = 0; i < nxOF; ++i)
            outStream << img[j * nxOF + i] << img[j * nxOF + i] << img[j * nxOF + i];
    outStream.close();
#endif

    // de-allocate optical flow solver
    delete[] img0;
    delete[] img1;

    delete[] ixix;
    delete[] ixiy;
    delete[] iyiy;

    delete[] u;
    delete[] v;
    delete[] uRHS;
    delete[] vRHS;
    delete[] uRes;
    delete[] vRes;

    delete[] uP;
    delete[] vP;
    delete[] uAP;
    delete[] vAP;

    // de-allocate ray tracer
    delete[] img;
    free(spheres);

    return 0;
}
