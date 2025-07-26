cmake_minimum_required(VERSION 3.15)
project(MatrixLib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include directory
include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(matrixlib STATIC
	src/matrix.cpp
)