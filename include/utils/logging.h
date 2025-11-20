#pragma once
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>

inline std::string current_time_str() {
  using namespace std::chrono;
  auto now = system_clock::now();
  auto t = system_clock::to_time_t(now);
  auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

  std::ostringstream oss;
  oss << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "."
      << std::setfill('0') << std::setw(3) << ms.count();
  return oss.str();
}

#define LOG_PRINT(level, stream_expr)                                    \
  do {                                                                   \
    std::cerr << "[" << level << "][" << current_time_str() << "] "      \
              << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") " \
              << stream_expr << std::endl;                               \
  } while (0)

#ifdef ENABLE_DEBUG_LOG
#define DEBUG(stream_expr) LOG_PRINT("DEBUG", stream_expr)
#else
#define DEBUG(stream_expr) \
  do {                     \
  } while (0)
#endif

#define INFO(stream_expr) LOG_PRINT("INFO", stream_expr)
#define WARN(stream_expr) LOG_PRINT("WARN", stream_expr)
#define ERROR(stream_expr) LOG_PRINT("ERROR", stream_expr)
