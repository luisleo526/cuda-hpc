cmake_minimum_required ( VERSION 3.25.2 )

project ( cuhpc LANGUAGES CXX CUDA )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )

find_package ( CUDAToolkit )
#message ( STATUS  "CUDAToolkit_FOUND=${CUDAToolkit_FOUND}" )
if ( CUDAToolkit_FOUND )
    enable_language( CUDA )
    set ( CMAKE_CUDA_STANDARD 20 )
    set ( CMAKE_CUDA_ARCHITECTURES "60;72;75" )
endif ()
message ( STATUS "CMAKE_CUDA_STANDARD=${CMAKE_CUDA_STANDARD}" )
find_package ( OpenMP REQUIRED )
#message ( STATUS  "OpenMP_FOUND=${OpenMP_FOUND}" )
list ( APPEND PRJ_LIBRARIES OpenMP::OpenMP_CXX )

set ( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=${OpenMP_CXX_FLAGS}" )

file(GLOB SOURCES "src/*.cu")
add_executable ( ${PROJECT_NAME} ${SOURCES} uccd_roc.cu)

# Make Tests

set_target_properties( ${PROJECT_NAME}
        PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        )

set_target_properties ( ${PROJECT_NAME}
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        )

target_include_directories ( ${PROJECT_NAME}
        PRIVATE
        ${PRJ_INCLUDE_DIRS}
        )

target_link_libraries( ${PROJECT_NAME}
        PRIVATE
        ${PRJ_LIBRARIES}
        )

target_compile_features ( ${PROJECT_NAME}
        PRIVATE
        ${PRJ_COMPILE_FEATURES}
        )