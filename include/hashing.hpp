#ifndef HASHING_IMPL_HPP
#define HASHING_IMPL_HPP

#include "../externals/ntHash/include/nthash/nthash.hpp"
#include "sequence_utils.h"
#include <ankerl/unordered_dense.h>
#include <iostream>
#include <string>

struct custom_key {
  const long pos; // trzeba template bo może być long
  const size_t hval;
};

class KmerHash {

  size_t k;

public:
  using is_avalanching = void;
  explicit KmerHash(size_t k) : k(k) {}
  size_t operator()(const custom_key key) const { return key.hval; }
};

class KmerEqual {
  size_t k;
  const char *sequence;

public:
  explicit KmerEqual(size_t k, const char *ptr) : k(k), sequence(ptr) {}

  bool operator()(const custom_key &ka, const custom_key &kb) const {
    std::string kmer_a;
    std::string kmer_b;
    if (ka.pos >= 0) {
      kmer_a = std::string_view(sequence + ka.pos, k);
    } else {
      kmer_a = get_reversed_strand(std::string_view(sequence + (-ka.pos), k));
    }

    if (kb.pos >= 0) {
      kmer_b = std::string_view(sequence + kb.pos, k);
    } else {
      kmer_b = get_reversed_strand(std::string_view(sequence + (-kb.pos), k));
    }
    return kmer_a == kmer_b;
  }
};

template <typename T>
void hash_sequences(const std::string &sequence, int k, T total_length,
                    std::vector<T> &out, T &c) {

  int shift = 0;
  const char *start = &sequence[1];
  const char *seq_ptr = sequence.data();
  nthash::BlindNtHash blind(start, 1, k, 0);

  auto kmer_hash = KmerHash(k);
  auto kmer_equal = KmerEqual(k, seq_ptr);
  auto kmers_dict =
      ankerl::unordered_dense::map<const custom_key, T, KmerHash, KmerEqual>(
          0, kmer_hash, kmer_equal);
  c = 1;
  size_t h;

  for (T i = 1; i < sequence.length() - k; i++) {
    const char *kmer_ptr = &sequence[i];
    h = blind.hashes()[0];
    T pos = i * (blind.get_forward_hash() >= blind.get_reverse_hash()) -
            i * (blind.get_forward_hash() < blind.get_reverse_hash());
    custom_key key = {pos, h};
    //  Skip kmers with '$' sign.
    if (sequence[i + k - 1] == '$') {
      shift = k - 1;
      blind.roll(sequence[i + k]);

      continue;
    }
    if (shift > 0) {
      shift--;
      blind.roll(sequence[i + k]);
      continue;
    }

    auto it = kmers_dict.find(key);
    if (pos > 0) {
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
  std::cout <<"Number of canonical k-mers: " <<kmers_dict.size() << "\n";
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
void get_reversed_index(const std::vector<T> &kmers_vec,
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
#endif // HASHING_IMPL_HPP
