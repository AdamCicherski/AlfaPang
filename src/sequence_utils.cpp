#include "sequence_utils.h"
#include <array>

void translate_sequence(std::string &sequence) {
  const std::unordered_map<char, char> translation_table{
      {'H', 'N'}, {'D', 'N'}, {'K', 'N'}, {'M', 'N'}, {'R', 'N'},
      {'S', 'N'}, {'W', 'N'}, {'Y', 'N'}, {'V', 'N'}, {'B', 'N'},
      {'a', 'A'}, {'c', 'C'}, {'t', 'T'}, {'g', 'G'}, {'n', 'N'},
      {'h', 'N'}, {'d', 'N'}, {'k', 'N'}, {'m', 'N'}, {'r', 'N'},
      {'s', 'N'}, {'w', 'N'}, {'y', 'N'}, {'v', 'N'}, {'b', 'N'}};
  for (char &c : sequence) {
    auto it = translation_table.find(c);
    if (it != translation_table.end()) {
      c = it->second;
    }
  }
}

std::string get_reversed_strand(const std::string_view &seq) {
  static char complement[256] = {0}; // Static array initialized only once

  if (complement['A'] == 0) {
    complement['A'] = 'T';
    complement['T'] = 'A';
    complement['C'] = 'G';
    complement['G'] = 'C';
    complement['N'] = 'N';
    complement['$'] = '$';
  }

  size_t n = seq.size();
  std::string rev_comp(n, ' ');
  size_t i = 0;
  for (; i < n; i++) {
    rev_comp[n - 1 - i] = complement[(unsigned char)seq[i]];
  }
  return rev_comp;
}

uint8_t base_hash(std::string s) {
  static constexpr std::array<uint8_t, 256> lookup = [] {
    std::array<uint8_t, 256> table = {};
    table['A'] = 1;
    table['C'] = 2;
    table['T'] = 3;
    table['G'] = 4;
    table['N'] = 5;
    table['$'] = 6;
    table['X'] = 7;
    return table;
  }();

  uint8_t i = lookup[s[0]];
  uint8_t j = lookup[s[1]];

  return i << 3 | j;
}


