if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -g -Wextra")
endif()


add_executable(
    Helmholtz Helmholtz.cc    
)

target_link_libraries(Helmholtz PRIVATE Matrix Grid Eigen Plot)
