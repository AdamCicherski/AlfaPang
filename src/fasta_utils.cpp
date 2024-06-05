#include "fasta_utils.h"
#include <fstream>

std::string read_sequences_from_fasta(const std::string &filename, std::vector<std::string> &names) {
    std::ifstream file(filename);
    std::string sequences = "$";
    std::string line, sequence;

    while (std::getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            auto name = line.substr(1);
            names.push_back(name);
        }
        if (line.empty() || line[0] == '>') {
            if (!sequence.empty()) {
                sequences += sequence + '$';
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }

    if (!sequence.empty()) {
        sequences += sequence + '$';
    }

    return sequences;
}
