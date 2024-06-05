cmake_minimum_required(VERSION 3.10)
project(AlfaPang)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add directories for includes and sources
include_directories(include)

# Add subdirectory for the source files
add_subdirectory(src)

# Add subdirectory for the external libraries
add_subdirectory(externals/emhash)
add_subdirectory(externals/robin-hood-hashing)

# Add executable
add_executable(AlfaPang main.cpp)

# Link the libraries
target_link_libraries(AlfaPang PRIVATE HashDefinitions FastaUtils SequenceUtils Hashing BFS GFAWriter emhash robin-hood-hashing)

