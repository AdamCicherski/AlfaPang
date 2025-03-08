#ifndef HASHING_IMPL_HPP
#define HASHING_IMPL_HPP

#include "../externals/emhash/hash_table8.hpp"
#include "../externals/ntHash/include/nthash/nthash.hpp"
#include "../externals/unordered_dense/include/ankerl/unordered_dense.h"
#include "sequence_utils.h"
#include <iostream>
#include <string>

struct custom_key {
  const char *ptr;
  const size_t hval;
};

class KmerHash2 {

  size_t k;

public:
  using is_avalanching = void;
  explicit KmerHash2(size_t k) : k(k) {}
  size_t operator()(const custom_key key) const { return key.hval; }
};

class KmerEqual2 {
  size_t k;

public:
  explicit KmerEqual2(size_t k) : k(k) {}

  bool operator()(custom_key ka, custom_key kb) const {
    // if (ka.hval != kb.hval) {
    //   return 0;
    // }
    std::string kmer_a(ka.ptr, k);
    std::string rev_comp_a = get_reversed_strand(kmer_a);
    std::string canonical_a = std::min(kmer_a, rev_comp_a);

    std::string kmer_b(kb.ptr, k);
    std::string rev_comp_b = get_reversed_strand(kmer_b);
    std::string canonical_b = std::min(kmer_b, rev_comp_b);

    return canonical_a == canonical_b;
  }
};

template <typename T>
void hash_sequences(const std::string &sequence, int k, T total_length,
                    std::vector<T> &out, T &c) {

  bool arr[256] = {0};
  arr[36] = 1;

  const unsigned num_hashes = 1;
  int shift = 0;
  const char *start = &sequence[1];
  nthash::BlindNtHash blind(start, 1, k, 0);
  auto kmer_hash = KmerHash2(k);
  auto kmer_equal = KmerEqual2(k);
  auto kmers_dict =
      ankerl::unordered_dense::map<const custom_key, T, KmerHash2, KmerEqual2>(
          0, kmer_hash, kmer_equal);
  c = 1;
  size_t h;

  for (T i = 1; i < sequence.length() - k; i++) {
    const char *kmer_ptr = &sequence[i];
    h = blind.hashes()[0];
    custom_key key = {kmer_ptr, h};
    std::string kmer(kmer_ptr, k);

    if (arr[sequence[i + k - 1]] == 1) {
      shift = k - 1;
      blind.roll(sequence[i + k]);

      continue;
    }
    if (shift > 0) {
      shift--;
      blind.roll(sequence[i + k]);

      continue;
    }

    const std::string reversed = get_reversed_strand(kmer);
    auto it = kmers_dict.find(key);
    if (reversed > kmer) {
      if (it != kmers_dict.end()) {
        out[i] = it->second;
      } else {
        kmers_dict[key] = c;
        out[i] = c;
        c++;
      }
    } else {
      if (it != kmers_dict.end()) {
        out[i] = -it->second;
      } else {
        kmers_dict[key] = c;
        out[i] = -c;
        c++;
      }
    }
    blind.roll(sequence[i + k]);
  }
  std::cout << kmers_dict.size() << "\n";
};


template <typename T>
void get_kmers_occ(const std::vector<T> &kmers_vec, std::vector<T> &kmers_occ) {
  for (size_t i = 0; i < kmers_vec.size(); i++) {

    kmers_occ[abs(kmers_vec[i])]++;
  }
  kmers_occ[0] = 0;
}

template <typename T> void cumulative_sum(std::vector<T> &vec) {
  for (size_t i = 1; i < vec.size(); i++) {
    vec[i] += vec[i - 1];
  }
}

template <typename T>
void get_kmers_pos(const std::vector<T> &kmers_vec,
                       const std::vector<T> &kmers_occ,
                       std::vector<T> &kmer_pos_map) {
  std::vector<T> used(kmers_occ.size(), 0);
  std::cout << kmer_pos_map.size() << "\n";
  for (size_t i = 0; i < kmers_vec.size(); i++) {
    size_t kmer = abs(kmers_vec[i]);
    if (kmer == 0) {
      continue;
    }
    kmer_pos_map[kmers_occ[kmer - 1] + used[kmer]] = i;
    used[kmer]++;
 }
}



void filter(const std::vector<std::vector<int>> &kmers_pos_map,
            std::vector<int> &kmers_vec, int max_occ) {
  for (auto positions : kmers_pos_map) {
    if (positions.size() > max_occ) {
      for (auto pos : positions) {
        kmers_vec[pos] = 0;
      }
    }
  }
}

#endif // HASHING_IMPL_HPP
