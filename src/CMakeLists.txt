add_library(HashDefinitions hash_definitions.cpp)
target_include_directories(HashDefinitions PRIVATE ../include)

add_library(FastaUtils fasta_utils.cpp)
target_include_directories(FastaUtils PRIVATE ../include)

add_library(SequenceUtils sequence_utils.cpp)
target_include_directories(SequenceUtils PRIVATE ../include)

add_library(Hashing INTERFACE)
target_include_directories(Hashing INTERFACE ../include)

add_library(BFS INTERFACE)
target_include_directories(BFS INTERFACE ../include)

add_library(GFAWriter INTERFACE)
target_include_directories(GFAWriter INTERFACE ../include)

add_library(GetStates INTERFACE)
target_include_directories(GetStates INTERFACE ../include)

add_library(ProcessBlock INTERFACE)
target_include_directories(ProcessBlock INTERFACE ../include)

add_library(CollapsePaths INTERFACE)
target_include_directories(CollapsePaths INTERFACE ../include)

