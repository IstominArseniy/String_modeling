cmake_minimum_required(VERSION 3.15)
project(String_modeling)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(String_modeling
        wav_generator.cpp)
