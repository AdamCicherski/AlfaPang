#ifndef GFA_WRITER_IMPL_HPP
#define GFA_WRITER_IMPL_HPP

#include "gfa_writer.hpp"
#include <fstream>

template <typename T1, typename T2>
void write_gfa(const std::vector<std::pair<T1, T1>> &labels,
               const std::vector<std::vector<T2>> &paths,
               const std::string &filename, const std::string &sequence,
               const std::vector<std::string> &names) {
  std::ofstream file(filename);
  for (size_t i = 0; i < paths.size(); ++i) {
    file << "S\t" << labels[i].first << "\t"
         << sequence.substr(labels[i].first,
                            labels[i].second - labels[i].first + 1)
         << "\n";
    for (size_t j = 0; j < paths[i].size(); ++j) {
      if (j > 0) {
        file << "L\t" << paths[i][j - 1] << "\t+\t" << paths[i][j]
             << "\t+\t0M\n";
      }
    }
  }
  file.close();
}

#endif // GFA_WRITER_IMPL_HPP
