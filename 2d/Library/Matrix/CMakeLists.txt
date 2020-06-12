#compile a library
include_directories(include)
file(GLOB SOURCES "src/*.cc")

if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -g -Wextra")
endif()


add_library(
    Matrix SHARED
    ${SOURCES}
)

target_include_directories(Matrix PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
