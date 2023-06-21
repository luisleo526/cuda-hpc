cmake_minimum_required(VERSION 3.0)

set(CMAKE_CUDA_ARCHITECTURES 75)

project(cuhpc CUDA)

set(CMAKE_CUDA_STANDARD 17)

file(GLOB SOURCES "src/*.h" "src/*.cpp" "src/*.cuh" "src/*.cu")

add_executable(cuhpc
        ${SOURCES}
        main.cu
        )

set_target_properties(cuhpc PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)

