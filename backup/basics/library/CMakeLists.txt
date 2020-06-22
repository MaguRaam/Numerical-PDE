cmake_minimum_required(VERSION 3.7)
set (CMAKE_CXX_STANDARD 11)

# Name of the project
project(weno)
include_directories(include)
file(GLOB SOURCES "source/*.cc")

# Add warnings 
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")
endif()

add_executable(a.out ${SOURCES}) 