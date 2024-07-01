#ifndef GFA_WRITER_IMPL_HPP
#define GFA_WRITER_IMPL_HPP

#include "gfa_writer.h"
#include <fstream>

template <typename T1, typename T2>
void write_gfa(const std::vector<std::pair<T1, T1>> &labels,
               const std::vector<std::vector<T2>> &paths,
               const std::string &filename, const std::string &sequence,
               const std::vector<std::string> &names) {
  robin_hood::unordered_set<robin_hood::pair<T2, T2>> edges;
  std::ofstream output;
  output.open(filename);
  output << "H\tVN:Z:1.0\n";
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
             << sign_right << "\t0M\n";
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

#endif // GFA_WRITER_IMPL_HPP
