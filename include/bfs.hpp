#ifndef BFS_IMPL_HPP
#define BFS_IMPL_HPP

#include <queue>
#include <string>
#include <vector>

// Norbert's trick + one vector
template <typename T1, typename T2>
void find_connected_bfs(const std::vector<T1> &kmers_vec,
                        const std::string &sequence,
                        const std::vector<T1> &reversed_index,
                        const std::vector<T1> &kmers_occ, int k,
                        std::vector<T2> &visited_nodes, T2 &v) {
  // needs a variable for bits shift depending on template
  int shift = sizeof(T1) * 8 - 1;
  for (size_t j = 0; j < sequence.size(); j++) {
    if (visited_nodes[j] != 0 || sequence[j] == '$') {
      continue;
    }
    std::queue<T1> queue;
    queue.push(j);
    visited_nodes[j] = v; // potrzebne dla filtra
    while (!queue.empty()) {
      auto pos = queue.front();
      queue.pop();
      for (int i = 0; i < k; i++) {
        if (pos - i >= 0) [[likely]] {
          auto kmer = kmers_vec[pos - i];
          if (kmer == 0) {
            continue;
          }

          auto kmer_index = abs(kmer);
          auto kmer_first_occ = reversed_index[kmers_occ[kmer_index - 1]];
          auto canonical_node =
              kmer_first_occ +
              ((kmer ^ kmers_vec[kmer_first_occ]) < 0) * (k - 1) +
              (1 | (kmer >> shift)) *
                  (1 | (kmers_vec[kmer_first_occ] >> shift)) * i;
          if (pos == canonical_node) {
            for (size_t occ_idx = kmers_occ[kmer_index - 1];
                 occ_idx < kmers_occ[kmer_index]; occ_idx++) {
              auto occ = reversed_index[occ_idx];

              auto new_pos = occ +
                             ((kmer ^ kmers_vec[occ]) < 0) * (k - 1) +
                             (1 | (kmer >> shift)) *
                                 (1 | (kmers_vec[occ] >> shift)) * i;
              if (visited_nodes[new_pos] != 0) {
                continue;
              } else {
                visited_nodes[new_pos] =
                    (sequence[new_pos] == sequence[j]) * v +
                    (sequence[new_pos] != sequence[j]) * (-v);
                queue.push(new_pos);
              }
            }
          } else {
            if (visited_nodes[canonical_node]) {
              continue;
            }
            queue.push(canonical_node);
            visited_nodes[canonical_node] =
                (sequence[canonical_node] == sequence[j]) * v +
                (sequence[canonical_node] != sequence[j]) * (-v);
          }
        } else [[unlikely]] {
          break;
        }
      }
    }
    v++;
  }
}

#endif // BFS_IMPL_HPP
