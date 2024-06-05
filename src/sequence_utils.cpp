#include "sequence_utils.h"
#include <algorithm>

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

std::string get_reversed_strand(const std::string_view &s) {
  static const char complement[] = {'T', 'G', 'A', 'C', 'N', '$'};
  std::string out;
  out.reserve(s.size());

  for (char i : s) {
    int index;
    switch (i) {
    case 'A':
      index = 0;
      break;
    case 'C':
      index = 1;
      break;
    case 'T':
      index = 2;
      break;
    case 'G':
      index = 3;
      break;
    case 'N':
      index = 4;
      break;
    case '$':
      index = 5;
      break;
    }
    out += complement[index];
  }

  std::reverse(out.begin(), out.end());
  return out;
}

char base_hash(std::string s) {
  char i;
  char j;
  switch (s[0]) {
  case 'A':
    i = 0;
    break;
  case 'C':
    i = 1;
    break;
  case 'T':
    i = 2;
    break;
  case 'G':
    i = 3;
    break;
  case 'N':
    i = 4;
    break;
  case '$':
    i = 5;
    break;
  case 'X':
    i = 6;
    break;
  }
  switch (s[1]) {
  case 'A':
    j = 0;
    break;
  case 'C':
    j = 1;
    break;
  case 'T':
    j = 2;
    break;
  case 'G':
    j = 3;
    break;
  case 'N':
    j = 4;
    break;
  case '$':
    j = 5;
    break;
  case 'X':
    j = 6;
    break;
  }
  return 7 * i + j + 1;
}
