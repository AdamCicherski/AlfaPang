#include "bfs.hpp"
#include "collapse_paths.hpp"
#include "fasta_utils.h"
#include "get_states.hpp"
#include "gfa_writer.hpp"
#include "hash_definitions.hpp"
#include "hashing.hpp"
#include "process_block.hpp"
#include "sequence_utils.h"
#include <climits>
#include <iostream>
#include <string>
#include <vector>
int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_fasta> <output_gfa> <k>\n";
    return 1;
  }

  std::string input = argv[1];
  std::string output_file = argv[2];
  int k = std::stoi(argv[3]);
  // Read .fa file
  std::cout << "Reading fasta" << std::endl;
  std::vector<std::string> names;
  std::string sequences = read_sequences_from_fasta(input, names);

  translate_sequence(sequences);

  long total_length = sequences.size();
  std::cout << "Total length of sequences: " << total_length - names.size() - 1
            << "\n";

  if (total_length < INT_MAX) {
    int total_length = sequences.size();
    std::cout << "Hashing k-mers..." << std::endl;
    std::vector<int> kmers_vec(total_length, 0);
    int kmers_number;
    hash_sequences(sequences, k, total_length, kmers_vec, kmers_number);

    std::vector<int> kmers_occ(kmers_number, 0);
    std::cout << "Counting k-mers occurrences" << std::endl;
    get_kmers_occ(kmers_vec, kmers_occ);
    cumulative_sum(kmers_occ);
    std::cout << "Building reversed index..." << std::endl;
    std::vector<int> kmers_pos(sequences.size(), 0);
    get_kmers_pos(kmers_vec, kmers_occ, kmers_pos);

    std::cout << "Starting BFS..." << std::endl;
    std::vector<int> choped(total_length, 0);
    int choped_nodes_number = 1;
    find_connected_bfs(kmers_vec, sequences, kmers_pos, kmers_occ, k, choped,
                       choped_nodes_number);
    kmers_vec = std::vector<int>();
    kmers_pos = std::vector<int>();

    std::cout << "Choped nodes number: " << choped_nodes_number << std::endl;
    std::cout << "Compacting unbranching paths..." << std::endl;

    std::vector<char> states(choped_nodes_number, 0);
    get_states(choped, sequences, states);
    std::vector<std::vector<int>> paths;
    std::vector<std::pair<int, int>> labels;
    collapse_paths<int, int>(states, choped, sequences, labels, paths);
    choped.clear();
    states.clear();
    std::cout << "Writing GFA..." << std::endl;

    write_gfa(labels, paths, output_file, sequences, names);

  } else {
    std::cout << "[Info] Using 64 bits ints for k-mers id\n";
    std::cout << "Hashing k-mers..." << std::endl;
    std::vector<long> kmers_vec(total_length, 0);
    long kmers_number;
    hash_sequences(sequences, k, total_length, kmers_vec, kmers_number);
    std::vector<long> kmers_occ(kmers_number, 0);
    std::cout << "Counting k-mers occurrences" << std::endl;
    get_kmers_occ(kmers_vec, kmers_occ);
    cumulative_sum(kmers_occ);
    std::cout << "Building reversed index..." << std::endl;
    std::vector<long> kmers_pos(sequences.size(), 0);
    get_kmers_pos(kmers_vec, kmers_occ, kmers_pos);

    if (kmers_number < INT_MAX) {
      int kmers_number = kmers_number;
      std::cout << "Starting BFS..." << std::endl;
      std::vector<int> choped(total_length, 0);
      int choped_nodes_number = 1;
      find_connected_bfs(kmers_vec, sequences, kmers_pos, kmers_occ, k, choped,
                         choped_nodes_number);
      kmers_vec = std::vector<long>();
      kmers_pos = std::vector<long>();

      std::cout << choped_nodes_number << std::endl;
      std::cout << "Compacting unbranched paths..." << std::endl;
      std::vector<char> states(choped_nodes_number, 0);
      get_states(choped, sequences, states);
      std::vector<std::vector<int>> paths;
      std::vector<std::pair<long, long>> labels;
      collapse_paths<long, int>(states, choped, sequences, labels, paths);
      choped.clear();
      states.clear();
      std::cout << "Writing GFA..." << std::endl;

      write_gfa(labels, paths, output_file, sequences, names);
    } else {
      std::vector<long> choped(total_length, 0);
      std::cout << "Starting BFS..." << std::endl;
      long choped_nodes_number = 1;
      find_connected_bfs(kmers_vec, sequences, kmers_pos, kmers_occ, k, choped,
                         choped_nodes_number);
      kmers_vec = std::vector<long>();
      kmers_pos = std::vector<long>();

      std::cout << choped_nodes_number << std::endl;
      std::cout << "Compacting unbranched paths..." << std::endl;
      std::vector<char> states(choped_nodes_number, 0);
      get_states(choped, sequences, states);
      std::vector<std::vector<long>> paths;
      std::vector<std::pair<long, long>> labels;
      collapse_paths<long, long>(states, choped, sequences, labels, paths);
      choped.clear();
      states.clear();
      std::cout << "Writing GFA..." << std::endl;

      write_gfa(labels, paths, output_file, sequences, names);
    }
  }
  return 0;
}
