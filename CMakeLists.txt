cmake_minimum_required(VERSION 3.16)

project(MinimaxComp)
set(CMAKE_CXX_COMPILER g++)

add_executable(degrees)
target_sources(degrees
	PRIVATE
		./run_compare.cpp
		./func.cpp
		./optimized_degrees.cpp
)
target_compile_options(degrees PRIVATE -Wall -O3 -std=c++11 -pthread -march=native)
target_link_libraries(degrees PRIVATE ntl gmp pthread m)
target_include_directories(degrees PUBLIC /usr/local/include /usr/local/include/NTL)
