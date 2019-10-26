cmake_minimum_required(VERSION 3.10.2)
project(CPPHomework1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -fsanitize=address --coverage")

add_executable(CPPHomework1 main.cpp Vector.cpp Matrix.cpp Base.h Vector.h Matrix.h catch.hpp)