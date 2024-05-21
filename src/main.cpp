#include <iostream>
#include <chrono>
#include <fmt/format.h>
#include <filesystem>
#include "Vec3.h"

int main()
{
    namespace fs=std::filesystem;
//    fs::create_directories("/transfer/particles/");
    std::cout<<"Particle System\n";
//    Emitter e(Vec3(0,0,0),10000);
//    for (int i=0; i<500; i++)
//    {
//        e.writeGeo(fmt::format("/transfer/particles/particle.{:04d}.geo", i));
//        e.update();
//    }
    return EXIT_SUCCESS;
}