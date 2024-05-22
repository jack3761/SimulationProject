#include <iostream>
#include <chrono>
#include <fmt/format.h>
#include <filesystem>
#include "Vec3.h"
#include "FluidSim.h"

int main()
{
    namespace fs=std::filesystem;
//    fs::create_directories("/transfer/particles/");
    std::cout<<"Particle System\n";
    FluidSim fluidSim(1000, 50, 50, 1, 1, 50);

    for (int i; i<100; i++)
    {
        float dt = 1/24;

        fluidSim.simulate(dt, -0.1, 0.9, 10, 1.9, true, true, 5.0, 5.0, 2.0, i);
    }
//    Emitter e(Vec3(0,0,0),10000);
//    for (int i=0; i<500; i++)
//    {
//        e.writeGeo(fmt::format("/transfer/particles/particle.{:04d}.geo", i));
//        e.update();
//    }
    return EXIT_SUCCESS;
}