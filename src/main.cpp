#include "bfs.hpp"
#include "fasta_utils.hpp"
#include "gfa_writer.hpp"
#include "hash_definitions.hpp"
#include "hashing.hpp"
#include "sequence_utils.hpp"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_fasta> <output_gfa> <k>\n";
    return 1;
  }

  std::string input_fasta = argv[1];
  std::string output_gfa = argv[2];
  int k = std::stoi(argv[3]);

  std::vector<std::string> names;
  std::string sequence = read_sequences_from_fasta(input_fasta, names);

  translate_sequence(sequence);

  using T = int; // Can be changed to int if required
  T total_length = static_cast<T>(sequence.size());
  std::vector<T> kmers_vec(total_length, 0);

  T c;
  hash_sequences(sequence, k, total_length, kmers_vec, c);

  std::vector<std::vector<T>> kmers_pos_map(c + 1);
  get_kmers_pos(kmers_vec, kmers_pos_map);

  std::vector<int> visited(total_length, 0);
  int v = 0;
  find_connected_bfs(kmers_vec, sequence, kmers_pos_map, k, visited, v);

  std::vector<std::pair<T, T>> nodes_labels;
  std::vector<std::vector<int>> paths;
  std::vector<char> states(total_length, 'N');

  get_states(kmers_vec, sequence, states);
  collapse_paths(states, visited, sequence, nodes_labels, paths);

  write_gfa(nodes_labels, paths, output_gfa, sequence, names);

  return 0;
}
