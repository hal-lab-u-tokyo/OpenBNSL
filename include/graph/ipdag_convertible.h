#pragma once
#include "graph/pdag.h"

/**
 * @brief Interface for converting internal graph structures to a unified PDAG.
 *
 * Structure learning algorithms may use different graph representations.
 * This interface ensures that all such graphs can be converted into a common
 * sparse PDAG form for downstream use (e.g., evaluation, comparison, or Python
 * interop).
 */
struct IPDAGConvertible {
  virtual ~IPDAGConvertible() = default;
  virtual PDAG to_pdag() const = 0;
};
