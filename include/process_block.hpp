#ifndef PROCESS_BLOCK_IMPL_HPP
#define PROCESS_BLOCK_IMPL_HPP

#include "../externals/robin-hood-hashing/src/include/robin_hood.h"
#include <string>
#include <vector>

template <typename T1, typename T2>
void process_block(
    robin_hood::unordered_node_map<robin_hood::pair<T2, T2>, T2> &nodes_id,
    std::vector<std::pair<T1, T1>> &nodes_labels, std::vector<T2> &path,
    T1 start, T1 stop, T2 &x, const std::string &sequence,
    const std::vector<T2> &transformed) {
  robin_hood::pair<T2, T2> n(transformed[start], transformed[stop]);
  robin_hood::pair<T2, T2> m(-n.second, -n.first);
  if (nodes_id.find(n) != nodes_id.end()) {
    path.push_back(nodes_id[n]);
  } else if (nodes_id.find(m) != nodes_id.end()) {
    path.push_back(-1 * nodes_id[m]);
  } else {
    nodes_id[n] = x;
    nodes_labels.push_back({start, stop + 1 - start});
    path.push_back(nodes_id[n]);
    x += 1;
  }
}

#endif // PROCESS_BLOCK_IMPL_HPP
