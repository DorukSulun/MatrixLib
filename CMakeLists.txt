cmake_minimum_required(VERSION 3.15)
project(MatrixLib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include directory
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/tests)	#for catch_amalgamated.hpp

add_library(matrixlib STATIC
	src/matrix.cpp
	src/utility.cpp
)


# TESTS
enable_testing()

add_executable(unit_tests
	tests/test_main.cpp
	tests/test_matrix_creation.cpp
	tests/test_utility.cpp
	tests/test_row_operations.cpp
	tests/test_arithmetic.cpp
	tests/test_transformations.cpp
	tests/test_special_matrices.cpp
	tests/test_pivot_tools.cpp
	tests/catch_amalgamated.cpp
)

target_include_directories(unit_tests PRIVATE tests)
target_link_libraries(unit_tests PRIVATE matrixlib)

add_test(NAME AllTests COMMAND unit_tests)