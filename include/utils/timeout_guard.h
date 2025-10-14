#pragma once
#include <atomic>
#include <chrono>
#include <iostream>

class TimeoutGuard {
 public:
  explicit TimeoutGuard(double timeout_sec = 0)
      : timeout_sec_(timeout_sec),
        start_(std::chrono::steady_clock::now()),
        timeout_flag_(false) {}

  bool expired() const {
    if (timeout_sec_ <= 0) return false;
    auto now = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(now - start_).count();
    if (elapsed > timeout_sec_) {
      timeout_flag_.store(true, std::memory_order_relaxed);
      return true;
    }
    return false;
  }

  bool is_timeout() const {
    return timeout_flag_.load(std::memory_order_relaxed);
  }

  double elapsed_sec() const {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::seconds>(now - start_)
        .count();
  }

 private:
  double timeout_sec_;
  std::chrono::steady_clock::time_point start_;
  mutable std::atomic<bool> timeout_flag_;
};
