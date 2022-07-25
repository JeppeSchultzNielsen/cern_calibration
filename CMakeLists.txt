cmake_minimum_required(VERSION 3.21)
project(calibration)

set(CMAKE_CXX_STANDARD 14)

find_package(AUSALIB REQUIRED)


add_executable(calibration main.cpp)
target_link_libraries(calibration PRIVATE AUSA::AUSA)
add_executable(efficiency efficiency.cpp)
target_link_libraries(efficiency PRIVATE AUSA::AUSA)
