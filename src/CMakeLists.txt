add_library(HashDefinitions hash_definitions.cpp)
target_include_directories(HashDefinitions PRIVATE ../include)

add_library(FastaUtils fasta_utils.cpp)
target_include_directories(FastaUtils PRIVATE ../include)

add_library(SequenceUtils sequence_utils.cpp)
target_include_directories(SequenceUtils PRIVATE ../include)

target_include_directories(Hashing PUBLIC ../include)
target_include_directories(BFS PUBLIC ../include)
target_include_directories(GFAWriter PUBLIC ../include)


