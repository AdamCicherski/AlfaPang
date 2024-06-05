#ifndef HASHING_HPP
#define HASHING_HPP

#include <vector>
#include <string>

template <typename T>
void hash_sequences(const std::string &sequence, int k, T total_length, std::vector<T> &out, T &c);

template <typename T>
void get_kmers_pos(const std::vector<T> &kmers_vec, std::vector<std::vector<T>> &kmers_pos_map);

#include "hashing_impl.hpp" // Include template implementations

#endif // HASHING_HPP

