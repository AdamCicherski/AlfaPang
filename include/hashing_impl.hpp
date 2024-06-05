#ifndef HASHING_IMPL_HPP
#define HASHING_IMPL_HPP

#include "../externals/emhash/hash_table8.hpp"
#include "sequence_utils.hpp"
#include <string>

template <typename T>
void hash_sequences(const std::string &sequence, int k, T total_length,
                    std::vector<T> &out, T &c) {
  emhash8::HashMap<std::string, T> kmers_dict;
  c = 1;
  for (T i = 0; i < sequence.length() - k; i++) {
    const std::string kmer = sequence.substr(i, k);
    if (kmer.find('$') < k) {
      continue;
    }
    const std::string reversed = get_reversed_strand(kmer);
    if (reversed > kmer) {
      auto it = kmers_dict.find(kmer);
      if (it != kmers_dict.end()) {
        out[i] = it->second;
      } else {
        kmers_dict[kmer] = c;
        out[i] = c;
        c++;
      }
    } else {
      auto it = kmers_dict.find(reversed);
      if (it != kmers_dict.end()) {
        out[i] = -it->second;
      } else {
        kmers_dict[reversed] = c;
        out[i] = -c;
        c++;
      }
    }
  }
}

template <typename T>
void get_kmers_pos(const std::vector<T> &kmers_vec,
                   std::vector<std::vector<T>> &kmers_pos_map) {
  for (T i = 0; i < kmers_vec.size(); i++) {
    if (kmers_vec[i] != 0) {
      kmers_pos_map[abs(kmers_vec[i])].push_back(i);
    }
  }
}

#endif // HASHING_IMPL_HPP
