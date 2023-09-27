#include <chrono>
#include <fstream>
#include <cstring>

#define USE_CIMG
#ifdef USE_CIMG

#include "Cimg.h"

#endif

#include "ray-tracing-util.h"
#include "intersect.h"
#include "optical-flow-base.h"
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
            // compute screen coordinates
            auto fov = 90;
            auto scale = tan(fov * 0.5 * M_PI / 180);
            auto imageAspectRatio = nx / (double) ny;

            auto origin = Vec3{0., 0., 2.};
//            auto origin = Vec3{2.1, 1.3, 1.7};

            auto px = (2 * (i + 0.5) / nx - 1) * scale * imageAspectRatio;
            auto py = (1 - 2 * (j + 0.5) / ny) * scale;
            auto direction = Vec3{px, py, -1};

            // rotate camera
            double alpha = t * 2. * M_PI, beta = 0., gamma = -0.25 * M_PI;

            origin = origin.rotate(alpha, beta, gamma);
            direction = direction.rotate(alpha, beta, gamma);

            // create and trace ray for current pixel
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

#ifdef USE_CIMG
    for (; t <= maxTime; t += dt) {
#else
        for (; t <= dt; t += dt) {
#endif
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

        for (size_t j = 0; j < nyOF; ++j)
            for (size_t i = 0; i < nxOF; ++i) {
                img0[j * nxOF + i] = 0;
                for (size_t jOff = 0; jOff < supersampling; ++jOff)
                    for (size_t iOff = 0; iOff < supersampling; ++iOff)
#ifdef USE_COLOR
                            for (auto dim = 0; dim < 3; ++dim)
                                img0[j * nxOF + i] += img[((supersampling * j + jOff) * nxRT + supersampling * i + iOff) * 3 + dim];
                img0[j * nxOF + i] /= 3;
#else
                img0[j * nxOF + i] += img[(supersampling * j + jOff) * nxRT + supersampling * i + iOff]
#endif
                img0[j * nxOF + i] /= supersampling * supersampling;
            }

        auto endMap = std::chrono::steady_clock::now();

        auto startOF = std::chrono::steady_clock::now();

        if (t > 0) {
            initGradientsAndRHS(nxOF, nyOF, img0, img1, mgIxIx[numLevels - 1], mgIxIy[numLevels - 1], mgIyIy[numLevels - 1], mgRhsU[numLevels - 1], mgRhsV[numLevels - 1]);
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
                    *imgShowOF.data(i, j, dim, 0) = std::sqrt(mgSolU[numLevels - 1][j * nxOF + i] * mgSolU[numLevels - 1][j * nxOF + i]
                                                              + mgSolV[numLevels - 1][j * nxOF + i] * mgSolV[numLevels - 1][j * nxOF + i]);

        imgShowRT.mirror('y');
        imgShowOF.mirror('y');
        imgShowRT.resize(nxOF - 2, nxOF - 2);
        imgShowOF.resize(nxOF - 2, nyOF - 2);
        static cimg_library::CImgDisplay displayRT(imgShowRT, "ray-tracing");
        displayRT = imgShowRT;
        static cimg_library::CImgDisplay displayOF(imgShowRT, "optical-flow");
        displayOF = imgShowOF;
        while (t >= maxTime && !displayRT.is_closed() && !displayOF.is_closed())
            displayRT.wait();
    }
#else
    }
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
