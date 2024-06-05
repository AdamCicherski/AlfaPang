#include "hash_definitions.h"

namespace robin_hood {

size_t hash<robin_hood::pair<int, int>>::operator()(
    robin_hood::pair<int, int> const &p) const noexcept {
  auto a = hash<int>{}(p.first);
  auto b = hash<int>{}(p.second);
  return hash_combine(a, b);
}

size_t hash<robin_hood::pair<int, int>>::hash_combine(size_t seed,
                                                      size_t v) noexcept {
  seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

size_t hash<robin_hood::pair<long, long>>::operator()(
    robin_hood::pair<long, long> const &p) const noexcept {
  auto a = hash<long>{}(p.first);
  auto b = hash<long>{}(p.second);
  return hash_combine(a, b);
}

size_t hash<robin_hood::pair<long, long>>::hash_combine(size_t seed,
                                                        size_t v) noexcept {
  seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

size_t hash<robin_hood::pair<long, int>>::operator()(
    robin_hood::pair<long, int> const &p) const noexcept {
  auto a = hash<long>{}(p.first);
  auto b = hash<int>{}(p.second);
  return hash_combine(a, b);
}

size_t hash<robin_hood::pair<long, int>>::hash_combine(size_t seed,
                                                       size_t v) noexcept {
  seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

} // namespace robin_hood
