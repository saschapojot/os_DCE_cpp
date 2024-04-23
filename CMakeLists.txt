cmake_minimum_required(VERSION 3.27)
project(os_DCE_cpp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("/home/polya/cppLocal/fftw/include")
include_directories("/home/polya/cppLocal/armadillo/include")
add_executable(os main.cpp
        evolution.cpp)

target_compile_options(os PRIVATE -march=native)


set_source_files_properties(
        evolution.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)
set_source_files_properties(
        main.cpp
        PROPERTIES
        COMPILE_OPTIONS "-O3"
)



set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/polya/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

find_package(Boost  REQUIRED filesystem json serialization)
if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()

target_link_libraries(os Boost::json)
target_link_libraries(os Boost::filesystem)
target_link_libraries(os Boost::serialization)

file(GLOB LIBS_SO "/home/polya/cppLocal/armadillo/lib/*.so")
target_link_libraries(os ${LIBS_SO})


file(GLOB LIBS_A "/home/polya/cppLocal/fftw/lib/*.a")
target_link_libraries(os ${LIBS_A})

