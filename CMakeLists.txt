cmake_minimum_required(VERSION 3.2)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

# name of the project It is best to use something different from the exe name
project(FluidSim)

# Here we set the C++ standard to use
set(CMAKE_CXX_STANDARD 17)

# Add include paths
include_directories(include)

# Now we add our target executable and the file it is built from.
add_executable(FluidSim
        include/Vec3.h
        include/FluidSim.h
        include/Particle.h
        include/FluidCell.h
        src/main.cpp
        src/FluidSim.cpp)

find_package(fmt CONFIG REQUIRED)
target_link_libraries(FluidSim PRIVATE fmt::fmt-header-only)
