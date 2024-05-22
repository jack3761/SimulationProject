#ifndef PARTICLE_H_
#define PARTICLE_H_
#include "Vec3.h"

struct Particle
{
    Particle() = default;
    Vec3 pos;
    Vec3 vel;
    Vec3 colour;
    float density;
    float restDensity;
};

#endif