#ifndef PROCESS_BLOCK_HPP
#define PROCESS_BLOCK_HPP

#include <vector>
#include <string>
#include "robin_hood.h"

template <typename T1, typename T2>
void process_block(
    robin_hood::unordered_node_map<robin_hood::pair<T2, T2>, T2> &nodes_id,
    std::vector<std::pair<T1, T1>> &nodes_labels, std::vector<T2> &path,
    T1 start, T1 stop, T2 &x, const std::string &sequence,
    const std::vector<T2> &transformed);

#include "process_block_impl.hpp"

#endif // PROCESS_BLOCK_HPP

