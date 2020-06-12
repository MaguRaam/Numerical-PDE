cmake_minimum_required(VERSION 3.7)
set (CMAKE_CXX_STANDARD 11)

# Name of the project
project(Diffusion)

# Add subdirectory:
add_subdirectory(Library) 

add_executable(
    Diffusion Diffusion.cc    
)

target_link_libraries(Diffusion PRIVATE Matrix Grid Eigen Plot) 
