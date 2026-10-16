#include "fasta_utils.h"
#include <fstream>
std::string read_sequences_from_fasta(const std::string &filename,
std::vector<std::string> &names) {
std::ifstream file(filename);
if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string sequences = "$";
    std::string line, seq_name, sequence;

    auto process_sequence = [&](const std::string &seq, const std::string
    &name) {
        size_t start = 0;
        while (start < seq.size()) {
            size_t n_start = seq.find('N', start);
            if (n_start == std::string::npos) {
                // brak N – cały fragment
                names.push_back(name + ":" + std::to_string(start) + "-" +
                std::to_string(seq.size())); sequences +=
                seq.substr(start) + '$'; break;
            }

            // fragment przed blokiem N
            if (n_start > start) {
                names.push_back(name + ":" + std::to_string(start) + "-" +
                std::to_string(n_start)); sequences += seq.substr(start,
                n_start - start) + '$';
            }

            // znajdź koniec bloku N
            size_t n_end = n_start;
            while (n_end < seq.size() && seq[n_end] == 'N')
                ++n_end;
            start = n_end;
        }
    };

    while (std::getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!sequence.empty()) {
                process_sequence(sequence, seq_name);
                sequence.clear();
            }
            seq_name = line.substr(1);
        } else if (!line.empty()) {
            sequence += line;
        }
    }

    if (!sequence.empty()) {
        process_sequence(sequence, seq_name);
    }

    return sequences;
}
