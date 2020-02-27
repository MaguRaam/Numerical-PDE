cmake_minimum_required(VERSION 3.7)
set (CMAKE_CXX_STANDARD 11)

# Name of the project
project(Helmholtz)
include_directories(include)
file(GLOB SOURCES "source/*.cc")

# Add warnings 
if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -g -Wextra")
endif()

add_executable(main ${SOURCES})
