#pragma once
#include <sstream>
#include <string>
#include <vector>

template <typename T>
std::string fmt_vec(const std::vector<T>& v) {
  std::ostringstream oss;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i > 0) oss << ",";
    oss << v[i];
  }
  return oss.str();
}