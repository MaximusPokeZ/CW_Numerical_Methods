cmake_minimum_required(VERSION 3.27)
project(CW_NM)

set(CMAKE_CXX_STANDARD 23)

include_directories(./include)
file(GLOB SRC src/*.cpp)

add_executable(CW_NM ${SRC})
