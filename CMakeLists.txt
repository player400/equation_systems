cmake_minimum_required(VERSION 3.27)
project(mn)

set(CMAKE_CXX_STANDARD 17)

add_executable(mn main.cpp
        Matrix.cpp
        Matrix.h
        gauss.h
        jacobi.h
        matrix_generators.h
        lu.h)
