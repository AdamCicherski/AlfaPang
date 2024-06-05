#ifndef BFS_IMPL_HPP
#define BFS_IMPL_HPP

#include "../externals/robin-hood-hashing/src/include/robin_hood.h"
#include "bfs.hpp"
#include <queue>
#include <vector>

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
}
#endif // BFS_IMPL_HPP
