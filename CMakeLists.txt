# forked from the template by j-levy
# How to use
# 1. mkdir build
# 2. cd build
# 3. cmake ..
# 4. make
# 5. <run your program>

cmake_minimum_required(VERSION 3.14)

set(PROJECT_NAME RAM-MSA)
set(DEBUG false)

project(${PROJECT_NAME} LANGUAGES CXX)

add_executable(${PROJECT_NAME} src/main.cpp
)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)

target_include_directories(${PROJECT_NAME} PRIVATE include)

if (DEBUG)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()
