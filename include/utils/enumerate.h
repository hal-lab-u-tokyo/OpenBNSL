#pragma once
#include <cstddef>
#include <iterator>
#include <utility>

namespace utils {

template <typename T>
auto enumerate(T&& iterable) {
  struct iterator {
    std::size_t i;
    decltype(std::begin(iterable)) iter;

    bool operator!=(const iterator& other) const { return iter != other.iter; }
    void operator++() {
      ++i;
      ++iter;
    }
    auto operator*() const { return std::pair{i, *iter}; }
  };

  struct wrapper {
    T& iterable;
    auto begin() { return iterator{0, std::begin(iterable)}; }
    auto end() { return iterator{0, std::end(iterable)}; }
  };

  return wrapper{iterable};
}

}  // namespace utils