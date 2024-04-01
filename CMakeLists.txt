cmake_minimum_required(VERSION 3.27)
project(os_DCE_cpp)

set(CMAKE_CXX_STANDARD 20)
include_directories("/home/polya/cppLocal/fftw/include")

add_executable(os main.cpp)

set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/polya/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)


find_package(Boost  REQUIRED filesystem json serialization)



file(GLOB LIBS_A "/home/polya/cppLocal/fftw/lib/*.a")


target_link_libraries(os ${LIBS_A})