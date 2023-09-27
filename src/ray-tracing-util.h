#pragma once

#include <iostream>
#include <cmath>

#include "vec3.h"

struct Ray {
    Vec3 origin;
    Vec3 direction;
};


#define USE_COLOR

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


void parseCLA(int argc, char *const *argv, size_t &numLevels, size_t &supersampling, size_t &mgIterations, double &dt, double &maxTime) {
    // default values
    numLevels = 10;
    supersampling = 2;
    mgIterations = 10;
    dt = 1e-2;
    maxTime = dt * 16;

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
