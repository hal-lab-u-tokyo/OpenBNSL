#pragma once
#include <iostream>

#define LOG_PRINT(level, stream_expr)                                        \
  do {                                                                       \
    std::cerr << "[" << level << "] " << __FILE__ << ":" << __LINE__ << " (" \
              << __func__ << ") " << stream_expr << std::endl;               \
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
