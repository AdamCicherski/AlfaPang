#ifndef FASTA_UTILS_HPP
#define FASTA_UTILS_HPP

#include <string>
#include <vector>

std::string read_sequences_from_fasta(const std::string &filename, std::vector<std::string> &names);

#endif // FASTA_UTILS_HPP
