#ifndef GFA_WRITER_HPP
#define GFA_WRITER_HPP

#include <string>
#include <vector>

template <typename T1, typename T2>
void collapse_paths(const std::vector<char> &states,
                    const std::vector<T2> &transformed,
                    const std::string &sequence,
                    std::vector<std::pair<T1, T1>> &nodes_labels,
                    std::vector<std::vector<T2>> &paths);

template <typename T1, typename T2>
void write_gfa(const std::vector<std::pair<T1, T1>> &labels,
               const std::vector<std::vector<T2>> &paths,
               const std::string &filename, const std::string &sequence,
               const std::vector<std::string> &names);

#include "gfa_writer_impl.hpp" // Include template implementations

#endif // GFA_WRITER_HPP
