cmake_minimum_required(VERSION 3.27)
project(libs C)

set(CMAKE_C_STANDARD 11)

add_library(libs STATIC ../../library.c
        vector.c
        vector.h
        vectorVoid.c
        vectorVoid.h)
