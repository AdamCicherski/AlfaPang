#ifndef SEQUENCE_UTILS_HPP
#define SEQUENCE_UTILS_HPP

#include <string>
#include <unordered_map>
#include <vector>

void translate_sequence(std::string &sequence);
std::string get_reversed_strand(const std::string_view &s);
char base_hash(std::string s);

#endif // SEQUENCE_UTILS_HPP

