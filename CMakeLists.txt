cmake_minimum_required(VERSION 3.10)
project(AlfaPang)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-march=native -mavx2 -include cstdint)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add project include directories
include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/externals/unordered_dense/include)

# Build ntHash directly with CMake
set(NTHASH_DIR ${CMAKE_SOURCE_DIR}/externals/ntHash)
include_directories(${NTHASH_DIR}/include)

file(GLOB NTHASH_SOURCES ${NTHASH_DIR}/src/*.cpp)

add_library(nthash STATIC ${NTHASH_SOURCES})
target_include_directories(nthash PUBLIC ${NTHASH_DIR}/include)

# Build other libraries
add_library(FastaUtils src/fasta_utils.cpp)
add_library(SequenceUtils src/sequence_utils.cpp)

# Build executable
add_executable(AlfaPang main.cpp)
target_link_libraries(AlfaPang PRIVATE FastaUtils SequenceUtils nthash)

