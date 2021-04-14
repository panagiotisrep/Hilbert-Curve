project(HilbertCurve)
cmake_minimum_required(VERSION 3.8)

set(CMAKE_BUILD_TYPE Debug)
add_definitions(${CMAKE_CXX_FLAGS} "-std=c++14")  # enable C++11 standard
add_definitions(${CMAKE_CXX_FLAGS} "-g")  # enable debuger

include_directories (BEFORE src)
include_directories(BEFORE include)

SET (SOURCE_FILES
        src/example.cpp
        )

SET (HEADER_FILES
        include/HilbertCurve.h)

add_executable (hilbertCurve ${SOURCE_FILES} ${HEADER_FILES})
