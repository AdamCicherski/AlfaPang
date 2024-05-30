
#include "emhash/hash_table8.hpp"
#include "robin-hood-hashing/src/include/robin_hood.h"
#include <iostream>
#include <vector>
// #include "absl/hash/hash.h"
#include <algorithm>
#include <climits>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
// #include "robin-map-1.2.1/include/tsl/robin_map.h"
#include <chrono> // For std::chrono::seconds
#include <deque>
#include <thread> // For std::this_thread::sleep_for
// #include "unordered_dense/src/ankerl.unordered_dense.cpp"

namespace robin_hood {

template <> struct hash<robin_hood::pair<int, int>> {
  size_t operator()(robin_hood::pair<int, int> const &p) const noexcept {
    auto a = hash<int>{}(p.first);
    auto b = hash<int>{}(p.second);
    return hash_combine(a, b);
  }

  static size_t hash_combine(size_t seed, size_t v) noexcept {
    // see
    // https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

template <> struct hash<robin_hood::pair<long, long>> {
  size_t operator()(robin_hood::pair<long, long> const &p) const noexcept {
    auto a = hash<long>{}(p.first);
    auto b = hash<long>{}(p.second);
    return hash_combine(a, b);
  }

  static size_t hash_combine(size_t seed, size_t v) noexcept {
    // see
    // https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

template <> struct hash<robin_hood::pair<long, int>> {
  size_t operator()(robin_hood::pair<long, int> const &p) const noexcept {
    auto a = hash<long>{}(p.first);
    auto b = hash<int>{}(p.second);
    return hash_combine(a, b);
  }

  static size_t hash_combine(size_t seed, size_t v) noexcept {
    // see
    // https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

} // namespace robin_hood

std::string read_sequences_from_fasta(const std::string &filename,
                                      std::vector<std::string> &names) {
  std::ifstream file(filename);
  std::string sequences = "$";
  std::string line, sequence;

  while (std::getline(file, line)) {

    if (!line.empty() && line[0] == '>') {
      auto name = line.substr(1);
      names.push_back(name);
    }
    if (line.empty() || line[0] == '>') {
      if (!sequence.empty()) {
        sequences += sequence + '$';
        sequence.clear();
      }
    } else {
      sequence += line;
    }
  }

  if (!sequence.empty()) {
    sequences += sequence + '$';
  }

  return sequences;
}

// Function to translate characters in a string based on a translation table
void translate_sequence(std::string &sequence) {

  const std::unordered_map<char, char> translation_table{
      {'H', 'N'}, {'D', 'N'}, {'K', 'N'}, {'M', 'N'}, {'R', 'N'},
      {'S', 'N'}, {'W', 'N'}, {'Y', 'N'}, {'a', 'A'}, {'c', 'C'},
      {'t', 'T'}, {'g', 'G'}, {'n', 'N'}};
  for (char &c : sequence) {
    auto it = translation_table.find(c);
    if (it != translation_table.end()) {
      c = it->second;
    }
  }
}

std::string get_reversed_strand(const std::string_view &s) {
  // Directly map characters to their complements using an array
  static const char complement[] = {'T', 'G', 'A', 'C', 'N', '$'};
  std::string out;
  out.reserve(s.size()); // Reserve space for the output string

  // Iterate over the input string and append the complement character to the
  // output string
  for (char i : s) {
    // Convert character to index ('A' -> 0, 'C' -> 1, etc.)
    int index;
    switch (i) {
    case 'A':
      index = 0;
      break;
    case 'C':
      index = 1;
      break;
    case 'T':
      index = 2;
      break;
    case 'G':
      index = 3;
      break;
    case 'N':
      index = 4;
      break;
    case '$':
      index = 5;
      break;
    }
    out += complement[index];
  }

  // Reverse the output string
  std::reverse(out.begin(), out.end());
  return out;
}

template <typename T>
void hash_sequences(const std::string &sequence, int k, T total_length,
                    std::vector<T> &out, T &c) {
  emhash8::HashMap<std::string, T> kmers_dict;
  c = 1;
  bool dbg = 0;
  for (T i = 0; i < sequence.length() - k; i++) {
    dbg = 0;
    const std::string kmer = sequence.substr(i, k);
    if (kmer.find('$') < k) {
      continue;
    }
    const std::string reversed = get_reversed_strand(kmer);
    if (dbg == 1) {
      std::cout << kmer << "\n" << reversed << "\n";
    }
    if (reversed > kmer) {
      auto it = (kmers_dict.find(kmer));
      if (it != kmers_dict.end()) {
        out[i] = it->second;
      } else {
        kmers_dict[kmer] = c;
        out[i] = c;
        c++;
      }
    } else {
      auto it = (kmers_dict.find(reversed));
      if (it != kmers_dict.end()) {
        out[i] = -it->second;
      } else {
        kmers_dict[reversed] = c;
        out[i] = -c;
        c++;
      }
    }
  }
  std::cout << "value of c : " << c << std::endl;
}

template <typename T>
void get_kmers_pos(const std::vector<T> &kmers_vec,
                   std::vector<std::vector<T>> &kmers_pos_map) {

  for (T i = 0; i < kmers_vec.size(); i++) {
    if (kmers_vec[i] != 0) {
      kmers_pos_map[abs(kmers_vec[i])].push_back(i);
    }
  }
}

// // Function to find connected components using BFS
template <typename T1, typename T2>
void find_connected_bfs(const std::vector<T1> &data,
                        const std::string &sequence,
                        const std::vector<std::vector<T1>> &kmers, int k,
                        std::vector<T2> &visited, T2 &v) {
  // needs a variable for bits shift depending on template
  for (T1 j = 0; j < sequence.size(); j++) {
    if (visited[j] == 0 && sequence[j] != '$') {
      auto seen_kmers = robin_hood::unordered_set<robin_hood::pair<T1, int>>{};
      std::vector<T1> queue = {j};
      while (!queue.empty()) {
        auto pos = queue.back();
        queue.pop_back();
        for (int i = 0; i < k; i++) {
          if (pos - i >= 0) [[likely]] {
            auto kmer_index = data[pos - i];
            int c = (kmer_index < 0) ? (k - 1 - i) : i;
            robin_hood::pair<T1, int> myVar(abs(kmer_index), c);
            if (seen_kmers.find(myVar) != seen_kmers.end()) {
              continue;
            }
            seen_kmers.insert(myVar);
            for (auto map : kmers[abs(kmer_index)]) {
              auto new_pos =
                  map + ((kmer_index ^ data[map]) < 0) * (k - 1) +
                  (1 | (kmer_index >> 31)) * (1 | (data[map] >> 31)) * i;
              if (visited[new_pos] != 0) {
                continue;
              } else {
                visited[new_pos] = (sequence[new_pos] == sequence[j]) ? v : -v;
                queue.push_back(new_pos);
              }
            }
          } else [[unlikely]] {
            break;
          }
        }
      }
      v++;
    }
  }
  std::cout << "The value of v is: " << v << std::endl;
}

char base_hash(std::string s) {
  char i;
  char j;
  switch (s[0]) {
  case 'A':
    i = 0;
    break;
  case 'C':
    i = 1;
    break;
  case 'T':
    i = 2;
    break;
  case 'G':
    i = 3;
    break;
  case 'N':
    i = 4;
    break;
  case '$':
    i = 5;
    break;
  case 'X':
    i = 6;
    break;
  }
  switch (s[1]) {
  case 'A':
    j = 0;
    break;
  case 'C':
    j = 1;
    break;
  case 'T':
    j = 2;
    break;
  case 'G':
    j = 3;
    break;
  case 'N':
    j = 4;
    break;
  case '$':
    j = 5;
    break;
  case 'X':
    j = 6;
    break;
  }
  return 7 * i + j + 1;
}

template <typename T>
void get_states(const std::vector<T> &choped, const std::string &sequence,
                std::vector<char> &states) {

  for (size_t i = 0; i < choped.size(); i++) {
    if (choped[i] == 0) {
      continue;
    }
    std::string s;
    s = std::string() + sequence[i - 1] + sequence[i + 1];
    if (choped[i] < 0) {
      s = get_reversed_strand(s);

      if (choped[i - 1] == -1 * choped[i]) {
        s[1] = 'X';
      }
      if (choped[i + 1] == -1 * choped[i]) {
        s[0] = 'X';
      }

    } else {

      if (choped[i - 1] == -1 * choped[i]) {
        s[0] = 'X';
      }
      if (choped[i + 1] == -1 * choped[i]) {
        s[1] = 'X';
      }
    }
    char state = base_hash(s);
    if (states[abs(choped[i])] == 0) // Not visited yet\n",
    {
      states[abs(choped[i])] = state;
    } else if (states[abs(choped[i])] == state) {
      continue;
    } else if (states[abs(choped[i])] == 49) {
      continue;
    }

    else if ((states[abs(choped[i])] - 1) / 7 != (state - 1) / 7 &&
             states[abs(choped[i])] % 7 != state % 7) {
      states[abs(choped[i])] = 49;
    } else if ((states[abs(choped[i])] - 1) / 7 != (state - 1) / 7) {
      states[abs(choped[i])] = base_hash(std::string() + 'X' + s[1]);
    } else if (states[abs(choped[i])] % 7 != state % 7) {
      states[abs(choped[i])] = base_hash(std::string() + s[0] + 'X');
    }
  }

  for (size_t i = 1; i < states.size();
       i++) {            // std::cout<<+states[i]<<std::endl;
    if (states[i] == 49) // Many inputs and Many outputs\n",
    {
      states[i] = 3;
      // std::cout<<i<<std::endl;
    } else if (states[i] % 7 == 0) // Many outputs\n",
    {
      // std::cout<<i<<std::endl;
      states[i] = 2;
      // std::cout<<2<<std::endl;
    } else if ((states[i] - 1) / 7 == 6) // Many inputs\n",
    {
      states[i] = 1;
      // std::cout<<i<<std::endl;
    } else // Non branching\n",
    {
      states[i] = 0;
    }
  }
}

template <typename T1, typename T2>
void process_block(
    robin_hood::unordered_node_map<robin_hood::pair<T2, T2>, T2> &nodes_id,
    std::vector<std::pair<T1, T1>> &nodes_labels, std::vector<T2> &path,
    T1 start, T1 stop, T2 &x, const std::string &sequence,
    const std::vector<T2> &transformed) { // std::cout<<transformed[start]<<"
                                          // "<< transformed[stop]<< std::endl;
  robin_hood::pair<T2, T2> n(transformed[start], transformed[stop]);
  robin_hood::pair<T2, T2> m(-n.second, -n.first);
  if (nodes_id.find(n) != nodes_id.end()) {
    path.push_back(nodes_id[n]);
  } else if (nodes_id.find(m) != nodes_id.end()) {
    path.push_back(-1 * nodes_id[m]);
  } else {
    nodes_id[n] = x;
    // std::string label = sequence[start:i+1-start];
    // nodes_labels.push_back(sequence.substr(start,stop+1-start));
    nodes_labels.push_back({start, stop + 1 - start});
    path.push_back(nodes_id[n]);
    x += 1;
  }
}

template <typename T1, typename T2>
void collapse_paths(const std::vector<char> &states,
                    const std::vector<T2> &transformed,
                    const std::string &sequence,
                    std::vector<std::pair<T1, T1>> &nodes_labels,
                    std::vector<std::vector<T2>> &paths) {

  robin_hood::unordered_node_map<robin_hood::pair<T2, T2>, T2> nodes_id;
  T1 start = 1;
  T2 x = 1;
  std::vector<T2> path;
  for (T1 i = 1; i < transformed.size(); i++) { // std::cout<<start<<std::endl;

    if ((states[abs(transformed[i])] == 1 && transformed[i] > 0) ||
        (states[abs(transformed[i])] == 2 &&
         transformed[i] < 0)) // BRANCH ON LEFT
    {
      if (i == start) {
        continue;
      }
      process_block(nodes_id, nodes_labels, path, start, i - 1, x, sequence,
                    transformed);
      start = i;

    }

    else if ((states[abs(transformed[i])] == 1 && transformed[i] < 0) ||
             (states[abs(transformed[i])] == 2 &&
              transformed[i] > 0)) // BRANCH ON RIGHT
    {                              // std::cout<<2<<std::endl;
      process_block(nodes_id, nodes_labels, path, start, i, x, sequence,
                    transformed);
      start = i + 1;
    } else if (states[abs(transformed[i])] == 3) // DOUBLE SIDE
    {                                            // std::cout<<3<<std::endl;
      if (i != start) {
        process_block(nodes_id, nodes_labels, path, start, i - 1, x, sequence,
                      transformed);
      }
      process_block(nodes_id, nodes_labels, path, i, i, x, sequence,
                    transformed);
      start = i + 1;
    }

    if (transformed[i] == 0) // ADD LAST BLOCK
    {
      if (start < i) {
        std::cout << "TO: " << i - 1 << std::endl;

        process_block(nodes_id, nodes_labels, path, start, i - 1, x, sequence,
                      transformed);
      }
      paths.push_back(path);
      path.clear();
      start = i + 1;
    }
  }
  std::cout << "Bloki: " << x << std::endl;
}

template <typename T1, typename T2>
void write_gfa(const std::vector<std::pair<T1, T1>> &labels,
               const std::vector<std::vector<T2>> &paths,
               const std::string &filename, const std::string &sequence,
               const std::vector<std::string> &names) {
  robin_hood::unordered_set<robin_hood::pair<T2, T2>> edges;
  std::ofstream output;
  output.open(filename);
  output << "H\n";
  for (size_t i = 0; i < labels.size(); i++) {
    output << "S\t" << i + 1 << "\t"
           << sequence.substr(labels[i].first, labels[i].second) << "\n";
  }

  for (auto path : paths) {
    for (size_t i = 0; i < path.size() - 1; i++) {
      robin_hood::pair<T2, T2> edge = {path[i], path[i + 1]};
      robin_hood::pair<T2, T2> reversed_edge = {-path[i + 1], -path[i]};
      if ((edges.find(edge) != edges.end()) ||
          (edges.find(reversed_edge) != edges.end())) {
        continue;
      }
      edges.insert(edge);
      char sign_left = (path[i] > 0) ? '+' : '-';
      char sign_right = (path[i + 1] > 0) ? '+' : '-';
      int left = abs(path[i]);
      int right = abs(path[i + 1]);

      output << "L\t" << left << "\t" << sign_left << "\t" << right << "\t"
             << sign_right << "\tOM\n";
    }
  }

  for (size_t i = 0; i < paths.size(); i++) {
    output << "P\t" << names[i] << "\t";

    for (size_t j = 0; j < paths[i].size() - 1; j++) {
      char sign = (paths[i][j] > 0) ? '+' : '-';
      output << abs(paths[i][j]) << sign << ',';
    }
    char sign = (paths[i][paths[i].size() - 1] > 0) ? '+' : '-';
    output << abs(paths[i][paths[i].size() - 1]) << sign << "\t*\n";
  }
}

int main(int argc, char *argv[]) {
  std::string input = argv[1];
  std::string output_file = argv[2];
  int k = std::stoi(argv[3]);


  // Read .fa file
  std::cout << "Czytam fasta" << std::endl;
  std::vector<std::string> names;
  std::string sequences = read_sequences_from_fasta(input, names);

  translate_sequence(sequences);

  long total_length = sequences.size();
  std::cout << "total len: " << total_length << "\n";

  if (total_length < INT_MAX) {
    int total_length = sequences.size();
    std::cout << "startuje hashowanie" << std::endl;
    std::vector<int> kmers_vec(total_length, 0);
    int kmers_number;
    hash_sequences(sequences, k, total_length, kmers_vec, kmers_number);
    std::cout << "mapuje kmery na sekwencje" << std::endl;
    std::vector<std::vector<int>> kmers_pos(kmers_number);
    get_kmers_pos(kmers_vec, kmers_pos);
    std::cout << "zaczynam bfs" << std::endl;
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
    std::cout << "licze states" << std::endl;
    std::vector<char> states(choped_nodes_number, 0);
    get_states(choped, sequences, states);
    std::vector<std::vector<int>> paths;
    std::vector<std::pair<int, int>> labels;
    collapse_paths<int, int>(states, choped, sequences, labels, paths);
    choped.clear();
    states.clear();
    write_gfa(labels, paths, output_file, sequences, names);

  } else {
    std::cout << "Potrzeba 8 bajtowych int\n";
    std::cout << "startuje hashowanie" << std::endl;
    std::vector<long> kmers_vec(total_length, 0);
    long kmers_number;
    hash_sequences(sequences, k, total_length, kmers_vec, kmers_number);
    std::cout << "mapuje kmery na sekwencje" << std::endl;
    std::vector<std::vector<long>> kmers_pos(kmers_number);
    get_kmers_pos(kmers_vec, kmers_pos);

    if (kmers_number < INT_MAX) {
      int kmers_number = kmers_number;
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
      std::cout << "licze states" << std::endl;
      std::vector<char> states(choped_nodes_number, 0);
      get_states(choped, sequences, states);
      std::vector<std::vector<int>> paths;
      std::vector<std::pair<long, long>> labels;
      collapse_paths<long, int>(states, choped, sequences, labels, paths);
      choped.clear();
      states.clear();
      write_gfa(labels, paths, output_file, sequences, names);
    } else {
      std::vector<long> choped(total_length, 0);
      long choped_nodes_number = 1;
      find_connected_bfs(kmers_vec, sequences, kmers_pos, k, choped,
                         choped_nodes_number);
      kmers_vec.clear();
      for (auto v : kmers_pos) {
        v.clear();
      }
      kmers_pos.clear();
      std::cout << choped_nodes_number << std::endl;
      std::cout << "licze states" << std::endl;
      std::vector<char> states(choped_nodes_number, 0);
      get_states(choped, sequences, states);
      std::vector<std::vector<long>> paths;
      std::vector<std::pair<long, long>> labels;
      collapse_paths<long, long>(states, choped, sequences, labels, paths);
      choped.clear();
      states.clear();
      write_gfa(labels, paths, output_file, sequences, names);
    }
  }
  return 0;
}
