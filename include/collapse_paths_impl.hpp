#ifndef COLLAPSE_PATHS_IMPL_HPP
#define COLLAPSE_PATHS_IMPL_HPP

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
  for (T1 i = 1; i < transformed.size(); i++) {

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
    {
      process_block(nodes_id, nodes_labels, path, start, i, x, sequence,
                    transformed);
      start = i + 1;
    } else if (states[abs(transformed[i])] == 3) // DOUBLE SIDE
    {
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

        process_block(nodes_id, nodes_labels, path, start, i - 1, x, sequence,
                      transformed);
      }
      paths.push_back(path);
      path.clear();
      start = i + 1;
    }
  }
}

#endif // COLLAPSE_PATHS_IMPL_HPP
