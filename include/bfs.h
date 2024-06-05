#ifndef BFS_HPP
#define BFS_HPP

#include <vector>
#include <string>

template <typename T1, typename T2>
void find_connected_bfs(const std::vector<T1> &data, const std::string &sequence, const std::vector<std::vector<T1>> &kmers, int k, std::vector<T2> &visited, T2 &v);

#include "bfs.hpp" // Include template implementations

#endif // BFS_HPP

