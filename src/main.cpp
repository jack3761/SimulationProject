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
    float dt = 1/60.0;
    int subSteps = 1;
    float flipRatio = 0.1;
    int numPressureIters = 50;
    int numParticleIters=2;

    // set up area to run sim
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
    int maxParticles = (numX * numY);

    FluidSim fluidSim(density, tankWidth, tankHeight, h, r, maxParticles);

    fluidSim.setNumParticles(maxParticles);

    float radius = std::min(tankWidth, tankHeight) / 2.0;
    float centerX = tankWidth /2.0;
    float centerY = tankHeight/2.0;
//    int numParticles = 1000;
    int p = 0;
    for (float rad = r; rad <= radius; rad += 2 * r) {
        int numPoints = static_cast<int>(2 * M_PI * rad / (2 * r));
        if (numPoints == 0) numPoints = 1; // Ensure there's at least one point for very small radii
        for (int n = 0; n < numPoints; ++n) {
            float angle = 2.0f * M_PI * n / numPoints;
            float xOffset = rad * cos(angle);
            float yOffset = rad * sin(angle);
            float particleX = centerX + xOffset;
            float particleY = centerY + yOffset;
            if (p < maxParticles) {
                fluidSim.setParticlePos(p * 2, particleX);
                fluidSim.setParticlePos(p * 2 + 1, particleY);
                ++p;
            } else {
                break;
            }
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

    int simIterations = 1500;
    for (int i=0; i<simIterations; i++)
    {
        if (i % (simIterations / 10) == 0) {
            int percentComplete = (i * 100) / simIterations;
            std::cout << percentComplete << "% complete\n";
        }
        fluidSim.simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters, 1.9, subSteps);
        fluidSim.writeGeo(fmt::format("/transfer/fluid/particle.{:04d}.geo", i));
    }

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = stop - start;

    std::cout << "Simulation complete in " << duration.count() << " seconds \n";
    return EXIT_SUCCESS;
}