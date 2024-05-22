#ifndef FLUIDCELL_H
#define FLUIDCELL_H

struct FluidCell
{
    float u; // Velocity component in x-direction
    float v; // Velocity component in y-direction
    float du; // Delta velocity component in x-direction
    float dv; // Delta velocity component in y-direction
    float prevU; // Previous velocity component in x-direction
    float prevV; // Previous velocity component in y-direction
    float p; // Pressure
    float s;
    int cellType;
};

#endif FLUIDCELL_H
