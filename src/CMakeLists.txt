# Create libraries from the source files

add_library(HashDefinitions hash_definitions.cpp)
target_include_directories(HashDefinitions PRIVATE ../include)

add_library(FastaUtils fasta_utils.cpp)
target_include_directories(FastaUtils PRIVATE ../include)

add_library(SequenceUtils sequence_utils.cpp)
target_include_directories(SequenceUtils PRIVATE ../include)

add_library(Hashing hashing.cpp)
target_include_directories(Hashing PRIVATE ../include)

add_library(BFS bfs.cpp)
target_include_directories(BFS PRIVATE ../include)

add_library(GFAWriter gfa_writer.cpp)
target_include_directories(GFAWriter PRIVATE ../include)

