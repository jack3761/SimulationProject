#include <iostream>
#include <chrono>
#include <fmt/format.h>
#include <filesystem>
#include "Vec3.h"
#include "FluidSim.h"

int main()
{
    namespace fs=std::filesystem;
    fs::create_directories("/transfer/fluid/");

    // set up scene

    float gravity = -9.81;
//    float dt = 1.0/120.0;
    float dt = 1/24.0;
    int numPressureIters = 50;
    int numParticleIters=2;

    float res = 100;

    float tankHeight = 10;
    float tankWidth = 15;
    float h = tankHeight / res;
    float density = 1000.0;

    float relWaterHeight = 0.8;
    float relWaterWidth = 0.6;

    float r = 0.3f * h;
    float dx = 2.0f * r;
    float dy = std::sqrt(3.0) / 2.0 * dx;

    int numX = static_cast<int>(std::floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx));
    int numY = static_cast<int>(std::floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy));
    int maxParticles = numX * numY;

    FluidSim fluidSim(density, tankWidth, tankHeight, h, r, maxParticles);

    fluidSim.setNumParticles(numX * numY);

//    fluidSim.setParticles();

    float p = 0;
    float tempPos = 0;
    for (size_t i=0; i<numX; i++)
    {
        for (size_t j=0; j<numY; j++)
        {
            tempPos = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
            fluidSim.setParticlePos(p++, tempPos);

            tempPos = h + r + dy * j;
            fluidSim.setParticlePos(p++, tempPos);
        }
    }

    int n = fluidSim.getfNumY();
    float s;
    for (size_t i=0; i<fluidSim.getfNumX(); i++)
    {
        for (size_t j=0; j<fluidSim.getfNumY(); j++)
        {
            s = 1.0;
            if (i==0 || i == fluidSim.getfNumX()-1 || j==0)
            {
                s = 0.0;
            }
            fluidSim.setS(i*n +j, s);
        }
    }

    std::cout << "Starting sim... \n";
    auto start = std::chrono::high_resolution_clock::now();

    int simIterations = 600;
    for (int i=0; i<simIterations; i++)
    {
        if (i % (simIterations / 10) == 0) {
            int percentComplete = (i * 100) / simIterations;
            std::cout << percentComplete << "% complete\n";
        }
        fluidSim.simulate(dt, -0.1, 0.9, 10, 1.9, true, true, 5.0, 5.0, 2.0, i);
        fluidSim.writeGeo(fmt::format("/transfer/fluid/particle.{:04d}.geo", i));
    }

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = stop - start;

    std::cout << "Simulation complete in " << duration.count() << " seconds \n";
    return EXIT_SUCCESS;
}