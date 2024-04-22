project ("DSM")

cmake_minimum_required(VERSION 3.5)

file(GLOB HEADER_FILES *.h)

file(GLOB SOURCE_FILES *.cpp)

add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_CXX_FLAGS} -std=c++17 -O3 -g -march=native -pthread -fopenmp")