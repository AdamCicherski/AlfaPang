#ifndef COLLAPSE_PATHS_HPP
#define COLLAPSE_PATHS_HPP

#include <vector>
#include <string>
#include "robin_hood.h"

template <typename T1, typename T2>
void collapse_paths(const std::vector<char> &states,
                    const std::vector<T2> &transformed,
                    const std::string &sequence,
                    std::vector<std::pair<T1, T1>> &nodes_labels,
                    std::vector<std::vector<T2>> &paths);

#include "collapse_paths_impl.hpp"

#endif // COLLAPSE_PATHS_HPP

