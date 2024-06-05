#ifndef BFS_IMPL_HPP
#define BFS_IMPL_HPP

#include "bfs.hpp"
#include <queue>
#include <vector>

template <typename T1, typename T2>
void find_connected_bfs(const std::vector<T1> &data,
                        const std::string &sequence,
                        const std::vector<std::vector<T1>> &kmers, int k,
                        std::vector<T2> &visited, T2 &v) {
  std::queue<size_t> q;
  for (size_t i = 0; i < data.size(); ++i) {
    if (!visited[i] && data[i] != 0) {
      q.push(i);
      visited[i] = ++v;
      while (!q.empty()) {
        size_t pos = q.front();
        q.pop();
        for (auto next : kmers[std::abs(data[pos])]) {
          if (!visited[next]) {
            visited[next] = v;
            q.push(next);
          }
        }
      }
    }
  }
}

#endif // BFS_IMPL_HPP
