#include "bfs.h"
#include "collapse_paths.h"
#include "fasta_utils.h"
#include "get_states.h"
#include "gfa_writer.h"
#include "hash_definitions.h"
#include "hashing.h"
#include "process_block.h"
#include "sequence_utils.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <climits>
#include <limits>
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
  std::cout << "Total length of sequences: " << total_length << "\n";

  if (total_length < INT_MAX) {
    int total_length = sequences.size();
    std::cout << "Hashing k-mers..." << std::endl;
    std::vector<int> kmers_vec(total_length, 0);
    int kmers_number;
    hash_sequences(sequences, k, total_length, kmers_vec, kmers_number);
    std::cout << "Building reversed index..." << std::endl;
    std::vector<std::vector<int>> kmers_pos(kmers_number);
    get_kmers_pos(kmers_vec, kmers_pos);
    std::cout << "Starting BFS..." << std::endl;
    std::vector<int> choped(total_length, 0);
    int choped_nodes_number = 1;
    find_connected_bfs(kmers_vec, sequences, kmers_pos, k, choped,
                       choped_nodes_number);
    kmers_vec.clear();
    for (auto v : kmers_pos) {
      v.clear();
    }
    kmers_pos.clear();
    std::cout << choped_nodes_number << std::endl;
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
    std::cout << "[Info] Using 64 bits ints for k-mers id";
    std::cout << "Hashing k-mers..." << std::endl;
    std::vector<long> kmers_vec(total_length, 0);
    long kmers_number;
    hash_sequences(sequences, k, total_length, kmers_vec, kmers_number);
    std::cout << "Building reversed index" << std::endl;
    std::vector<std::vector<long>> kmers_pos(kmers_number);
    get_kmers_pos(kmers_vec, kmers_pos);

    if (kmers_number < INT_MAX) {
      int kmers_number = kmers_number;
      std::vector<int> choped(total_length, 0);
      int choped_nodes_number = 1;
      std::cout << "Starting BFS..." << std::endl;
      find_connected_bfs(kmers_vec, sequences, kmers_pos, k, choped,
                         choped_nodes_number);
      kmers_vec.clear();
      for (auto v : kmers_pos) {
        v.clear();
      }
      kmers_pos.clear();
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
      long choped_nodes_number = 1;
      std::cout << "Starting BFS..." << std::endl;
      find_connected_bfs(kmers_vec, sequences, kmers_pos, k, choped,
                         choped_nodes_number);
      kmers_vec.clear();
      for (auto v : kmers_pos) {
        v.clear();
      }
      kmers_pos.clear();
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
