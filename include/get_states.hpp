#ifndef GET_STATES_IMPL_HPP
#define GET_STATES_IMPL_HPP

#include "sequence_utils.h"
#include <string>
#include <vector>

template <typename T>
void get_states(const std::vector<T> &choped, const std::string &sequence,
                std::vector<uint8_t> &states) {

  for (size_t i = 0; i < choped.size(); i++) {
    if (choped[i] == 0) {
      continue;
    }
    std::string s;
    s = std::string() + sequence[i - 1] + sequence[i + 1];
    if (choped[i] < 0) {
      s = get_reversed_strand(s);

      if (choped[i - 1] == -1 * choped[i]) {
        s[1] = 'X';
      }
      if (choped[i + 1] == -1 * choped[i]) {
        s[0] = 'X';
      }

    } else {

      if (choped[i - 1] == -1 * choped[i]) {
        s[0] = 'X';
      }
      if (choped[i + 1] == -1 * choped[i]) {
        s[1] = 'X';
      }
    }
    uint8_t state = base_hash(s);
    uint8_t new_l = state >> 3;
    uint8_t new_r = state & 7;
    auto node_id = abs(choped[i]);
    uint8_t old_l = states[node_id] >> 3;
    uint8_t old_r = states[node_id] & 7;

    if (states[node_id] == 0) // Not visited yet\n",
    {
      states[node_id] = state;
    } else if (states[node_id] == state) {
      continue;
    } else if (states[node_id] == 63) {
      continue;
    } else if (old_l != new_l && old_r != new_r) {
      states[node_id] = 63;
    } else if (old_l != new_l) {
      states[node_id] = base_hash(std::string() + 'X' + s[1]);
    } else if (old_r != new_r) {
      states[node_id] = base_hash(std::string() + s[0] + 'X');
    }
  }
 constexpr uint8_t LUT[256] = {
        // 
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2,
        0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 3,
        // Remaining values default to 0
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
  for (size_t i = 1; i < states.size(); i++) {
     states[i]=LUT[states[i]];}
}


#endif // GET_STATES_IMPL_HPP
