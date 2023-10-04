#pragma once

#include <iostream>
#include <cmath>

#include "vec3.h"


struct Ray {
    Vec3 origin;
    Vec3 direction;
};


//#define USE_COLOR

#ifdef USE_COLOR
using Color = Vec3;
#else
using Color = double;
#endif


struct Sphere {
    Vec3 pos;
    double rSq;     // radius squared

    Color diffuseColor;
    Color specularColor;
    double specularExp;
    Color reflectionColor;
};

//Sphere spheres[] = {
//        {{0.,        0.,        -160.},    128 * 128.,            {0.8}, {0.2}, 16, {0.8}},
//        {{0.,        0.,        0.},       0.4 * 0.4,             {0.8}, {0.2}, 8,  {0.2}},
//        {{0.272166,  0.272166,  0.544331}, (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{0.643951,  0.172546,  0.},       (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{0.172546,  0.643951,  0.},       (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{-0.371785, 0.099620,  0.544331}, (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{-0.471405, 0.471405,  0.},       (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{-0.643951, -0.172546, 0.},       (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{0.099620,  -0.371785, 0.544331}, (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{-0.172546, -0.643951, 0.},       (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}},
//        {{0.471405,  -0.471405, 0.},       (1. / 6.) * (1. / 6.), {0.2}, {0.8}, 32, {0.8}}
//};
//
//constexpr auto numSpheres = sizeof(spheres) / sizeof(Sphere);

Sphere *spheres;
int numSpheres;


inline void initSpheres(double t) {
    spheres[0] = {{0., 0., -130.}, 128 * 128., Color{0.8}, Color{0.2}, 16, Color{0.8}};
    spheres[1] = {{0., 0., 0.}, 0.5 * 0.5, Color{0.8}, Color{0.2}, 8, Color{0.2}};

    for (auto i = 2; i < numSpheres; ++i) {
        auto r = 1. / 8.;

        auto alpha = 2. * M_PI * (-(t + (double) (i - 2) / (numSpheres - 2)));
        auto beta = 2. * M_PI * ((i - 2 >= (numSpheres - 2) / 2 ? -1 : 1) / 16.);
        auto gamma = 2. * M_PI * (0.);

        auto pos = Vec3{sqrt(spheres[1].rSq) + ((numSpheres - 2) / 4.) * r, 0, 0}.rotate(alpha, beta, gamma);

        spheres[i] = {pos, r * r, Color{0.2}, Color{0.8}, 32, Color{0.8}};
    }
}


struct Light {
    Vec3 position; // 3D position
    Color color;
};
//using Light = Vec3;

#ifdef USE_COLOR
Light lights[] = {
//        {{4.,  3.,  2.}, {1., 0., 0.}},
//        {{1.,  -4., 4.}, {0., 1., 0.}},
//        {{-3., 1.,  5.}, {0., 0., 1.}}
        {{-16., -16., 16.}, {1., 0., 0.}},
        {{16.,  -16., 16.}, {0., 1., 0.}},
        {{0.,   16.,  16.}, {0., 0., 1.}}
//        {{-16., -16., 16.}, {1., 1., 1.}},
//        {{16.,  -16., 16.}, {1., 1., 1.}},
//        {{0.,   16.,  16.}, {1., 1., 1.}}
};
#else
Light lights[] = {
//        {{4.,  3.,  2.}, 1./3.},
//        {{1.,  -4., 4.}, 1./3.},
//        {{-3., 1.,  5.}, 1./3.}
        {{-16., -16., 16.}, 1. / 3.},
        {{16.,  -16., 16.}, 1. / 3.},
        {{0.,   16.,  16.}, 1. / 3.}
};
#endif

constexpr auto numLights = sizeof(lights) / sizeof(Light);
