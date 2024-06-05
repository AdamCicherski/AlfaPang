#ifndef GET_STATES_HPP
#define GET_STATES_HPP

#include <vector>
#include <string>
#include "sequence_utils.hpp"

template <typename T>
void get_states(const std::vector<T> &choped, const std::string &sequence,
                std::vector<char> &states);

#include "get_states_impl.hpp"

#endif // GET_STATES_HPP

