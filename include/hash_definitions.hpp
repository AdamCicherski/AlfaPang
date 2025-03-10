#ifndef HASH_DEFINITIONS_HPP
#define HASH_DEFINITIONS_HPP

#include "../externals/robin-hood-hashing/src/include/robin_hood.h"

namespace robin_hood {

template <typename T1, typename T2>
struct hash<robin_hood::pair<T1, T2>> {
  size_t operator()(robin_hood::pair<T1, T2> const &p) const noexcept {
    return hash_combine(std::hash<T1>{}(p.first), std::hash<T2>{}(p.second));
  }

  static size_t hash_combine(size_t seed, size_t v) noexcept {
    return seed ^ (v + 0x9e3779b9 + (seed << 6) + (seed >> 2));
  }
};

} // namespace robin_hood

#endif // HASH_DEFINITIONS_HPP

