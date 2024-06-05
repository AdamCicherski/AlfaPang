#ifndef HASH_DEFINITIONS_HPP
#define HASH_DEFINITIONS_HPP

#include "externals/robin-hood-hashing/src/include/robin_hood.h"

namespace robin_hood {

template <> struct hash<robin_hood::pair<int, int>> {
  size_t operator()(robin_hood::pair<int, int> const &p) const noexcept;
  static size_t hash_combine(size_t seed, size_t v) noexcept;
};

template <> struct hash<robin_hood::pair<long, long>> {
  size_t operator()(robin_hood::pair<long, long> const &p) const noexcept;
  static size_t hash_combine(size_t seed, size_t v) noexcept;
};

template <> struct hash<robin_hood::pair<long, int>> {
  size_t operator()(robin_hood::pair<long, int> const &p) const noexcept;
  static size_t hash_combine(size_t seed, size_t v) noexcept;
};

} // namespace robin_hood

#endif // HASH_DEFINITIONS_HPP
