#ifndef GET_STATES_IMPL_HPP
#define GET_STATES_IMPL_HPP

template <typename T>
void get_states(const std::vector<T> &choped, const std::string &sequence,
                std::vector<char> &states) {

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
    char state = base_hash(s);
    if (states[abs(choped[i])] == 0) // Not visited yet\n",
    {
      states[abs(choped[i])] = state;
    } else if (states[abs(choped[i])] == state) {
      continue;
    } else if (states[abs(choped[i])] == 49) {
      continue;
    }

    else if ((states[abs(choped[i])] - 1) / 7 != (state - 1) / 7 &&
             states[abs(choped[i])] % 7 != state % 7) {
      states[abs(choped[i])] = 49;
    } else if ((states[abs(choped[i])] - 1) / 7 != (state - 1) / 7) {
      states[abs(choped[i])] = base_hash(std::string() + 'X' + s[1]);
    } else if (states[abs(choped[i])] % 7 != state % 7) {
      states[abs(choped[i])] = base_hash(std::string() + s[0] + 'X');
    }
  }

  for (size_t i = 1; i < states.size(); i++) {
    if (states[i] == 49) // Many in and many out",
    {
      states[i] = 3;
      // std::cout<<i<<std::endl;
    } else if (states[i] % 7 == 0) // Many out",
    {
      // std::cout<<i<<std::endl;
      states[i] = 2;
      // std::cout<<2<<std::endl;
    } else if ((states[i] - 1) / 7 == 6) // Many in",
    {
      states[i] = 1;
    } else // Non branching",
    {
      states[i] = 0;
    }
  }
}

#endif // GET_STATES_IMPL_HPP
