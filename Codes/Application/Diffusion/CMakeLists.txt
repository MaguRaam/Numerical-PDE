if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -g -Wextra")
endif()


add_executable(
    Diffusion Diffusion.cc    
)

target_link_libraries(Diffusion PRIVATE Matrix Grid Eigen Plot)
