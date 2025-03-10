#ifndef BFS_IMPL_HPP
#define BFS_IMPL_HPP

#include <string>
#include <queue>
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
  for (T1 j = 0; j < sequence.size(); j++) {
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
          auto c = reversed_index[kmers_occ[kmer_index - 1]];
          auto new_edge =
              c + ((kmer ^ kmers_vec[c]) < 0) * (k - 1) +
              (1 | (kmer >> shift)) * (1 | (kmers_vec[c] >> shift)) * i;
          if (pos == new_edge) {
            for (size_t occ = kmers_occ[kmer_index - 1];
                 occ < kmers_occ[kmer_index]; occ++) {
              auto map = reversed_index[occ];

              auto new_pos =
                  map + ((kmer ^ kmers_vec[map]) < 0) * (k - 1) +
                  (1 | (kmer >> shift)) * (1 | (kmers_vec[map] >> shift)) * i;
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
            if (visited_nodes[new_edge]) {
              continue;
            }
            queue.push(new_edge);
            visited_nodes[new_edge] =
                (sequence[new_edge] == sequence[j]) * v +
                (sequence[new_edge] != sequence[j]) * (-v);
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
