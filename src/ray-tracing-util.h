#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>

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
    Vec3 initPos;   // initial position before moving
    Vec3 pos;       // current position
    double rSq;     // radius squared

    Color diffuseColor;
    Color specularColor;
    double specularExp;
    Color reflectionColor;
};


Sphere *spheres;
int numSpheres;

constexpr auto boundaryExtentXZ = 2.; // extent in each direction (positive & negative) in the x- and z-axes
constexpr auto boundaryExtentY = 2 * boundaryExtentXZ;


inline void initSpheres() {
    numSpheres = 48;
    spheres = static_cast<Sphere *>(malloc(numSpheres * sizeof(Sphere)));
    std::mt19937 randGen(0);// constant seed
    std::uniform_real_distribution<> posXZDistribution(-boundaryExtentXZ, boundaryExtentXZ);
    std::uniform_real_distribution<> posYDistribution(-boundaryExtentY, boundaryExtentY);
    std::uniform_real_distribution<> rDistribution(0.3, 0.6);
    for (auto s = 0; s < numSpheres; ++s) {
        auto pos = Vec3{posXZDistribution(randGen), posYDistribution(randGen), posXZDistribution(randGen)};
        auto r = rDistribution(randGen);

        auto invalidPos = false;

        // check for side walls
        if (std::sqrt(pos.x * pos.x + pos.z * pos.z) + r > boundaryExtentXZ * (0.625 + 0.375 * sin((pos.y / boundaryExtentY + 0.5) * M_PI)))
            invalidPos = true;

        // check for collisions with other particles; also check periodic mirror particle positions
        for (auto sCmp = 0; sCmp < s && !invalidPos; ++sCmp)
            if (std::min({
                                 (spheres[sCmp].pos - pos).magnitude(),
                                 (spheres[sCmp].pos - pos + Vec3{0., 2 * boundaryExtentY, 0.}).magnitude(),
                                 (spheres[sCmp].pos - pos + Vec3{0., -2 * boundaryExtentY, 0.}).magnitude()
                         }) < 1.25 * (std::sqrt(spheres[sCmp].rSq) + r))
                invalidPos = true;

        if (invalidPos) {
            --s;
            continue;
        }

        spheres[s] = {pos, pos, r * r, Color{0.5}, Color{0.8}, 32, Color{0.8}};
    }
}

inline void updateSpherePositions(double t) {
    for (auto s = 0; s < numSpheres; ++s) {
        spheres[s].pos.y = spheres[s].initPos.y - t * 2 * boundaryExtentY;
        while (spheres[s].pos.y < -boundaryExtentY)
            spheres[s].pos.y += 2 * boundaryExtentY;
    }
}


struct Light {
    Vec3 position; // 3D position
    Color color;
};
//using Light = Vec3;

#ifdef USE_COLOR
Light lights[] = {
        {{-16., 0., -16.}, {1.,  0.5, 0.5}},
        {{16.,  0., -16.}, {0.5, 1.,  0.5}},
        {{0.,   0., 16.},  {0.5, 0.5, 1.}}
};
#else
Light lights[] = {
        {{-16., 0., -16.}, 1. / 3.},
        {{16.,  0., -16.}, 1. / 3.},
        {{0.,   0., 16.},  1. / 3.}
};
#endif

constexpr auto numLights = sizeof(lights) / sizeof(Light);
