cmake_minimum_required(VERSION 3.16)
project(sim-particule)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -g -O3 -march=native")

set(CMAKE_CXX_STANDARD 17)

set(MAIN_SOURCE_FILES
    src/main.cpp
    src/part.cpp
    src/utils.cpp
)

set(HEADER_FILES
    include/part.hpp
    include/utils.hpp
)

add_executable(spart ${MAIN_SOURCE_FILES})

target_include_directories(spart PRIVATE ${PROJECT_SOURCE_DIR}/include)