include_directories(include)
file(GLOB SOURCES "src/*.cc")

if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -g -Wextra")
endif()

add_library(
    Eigen SHARED
    ${SOURCES}
)

target_link_libraries(Eigen PRIVATE Matrix Grid)
target_include_directories(Eigen PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
