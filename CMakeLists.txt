cmake_minimum_required(VERSION 3.10)
project(AlfaPang)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-march=native -mavx2)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/externals/unordered_dense/include)
# Build ntHash
set(NTHASH_DIR ${CMAKE_SOURCE_DIR}/externals/ntHash)
set(NTHASH_BUILD_DIR ${NTHASH_DIR}/build)
set(NTHASH_INSTALL_DIR ${NTHASH_DIR}/install)

# command to build ntHash using Meson
add_custom_target(build_nthash
    COMMAND meson setup --buildtype=release --prefix=${NTHASH_INSTALL_DIR} ${NTHASH_BUILD_DIR}
    COMMAND meson install -C ${NTHASH_BUILD_DIR}
    WORKING_DIRECTORY ${NTHASH_DIR}
    COMMENT "Building ntHash with Meson..."
)

include_directories(${NTHASH_INSTALL_DIR}/include)
link_directories(${NTHASH_INSTALL_DIR}/lib)
add_library(FastaUtils src/fasta_utils.cpp)
add_library(SequenceUtils src/sequence_utils.cpp)



add_executable(AlfaPang main.cpp)
target_link_libraries(AlfaPang PRIVATE FastaUtils SequenceUtils nthash)


add_dependencies(AlfaPang build_nthash)

