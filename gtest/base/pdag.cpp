#include "base/PDAG.h"

#include <gtest/gtest.h>

#include <algorithm>

// Test case 1: Default Construction - verify that no edge exists initially.
TEST(PDAGTest, DefaultConstruction) {
  constexpr size_t n = 5;
  PDAG g(n);
  // Initially, every edge should be absent.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      EXPECT_FALSE(g.has_edge(i, j))
          << "Edge unexpectedly exists from " << i << " to " << j;
    }
  }
}

// Test case 2: Adding and removing edges along with boundary checking.
TEST(PDAGTest, AddRemoveEdge) {
  constexpr size_t n = 5;
  PDAG g(n);

  // Add an edge from node 1 to node 3.
  g.add_edge(1, 3);
  EXPECT_TRUE(g.has_edge(1, 3)) << "Edge from 1 to 3 was not added.";

  // Edge from 3 to 1 does not exist unless explicitly added.
  EXPECT_FALSE(g.has_edge(3, 1)) << "Edge from 3 to 1 should not exist.";

  // Adding the same edge again should throw an exception.
  EXPECT_THROW(g.add_edge(1, 3), std::invalid_argument);

  // Remove the edge and confirm its removal.
  g.remove_edge(1, 3);
  EXPECT_FALSE(g.has_edge(1, 3)) << "Edge from 1 to 3 was not removed.";

  // Removing a non-existent edge should throw an exception.
  EXPECT_THROW(g.remove_edge(1, 3), std::invalid_argument);

  // Out-of-range index accesses should throw exceptions.
  EXPECT_THROW(g.has_edge(5, 0), std::out_of_range);
  EXPECT_THROW(g.add_edge(0, 5), std::out_of_range);
  EXPECT_THROW(g.remove_edge(5, 5), std::out_of_range);
}

// Test case 3: complete_graph - should create a complete graph without
// self-loops.
TEST(PDAGTest, CompleteGraph) {
  constexpr size_t n = 5;
  PDAG g(n);
  g.complete_graph();

  // Self-loops should be removed.
  for (size_t i = 0; i < n; ++i) {
    EXPECT_FALSE(g.has_edge(i, i))
        << "Self-loop at node " << i << " was not removed.";
  }
  // Every distinct pair (i, j) where i != j should have an edge.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i != j) {
        EXPECT_TRUE(g.has_edge(i, j))
            << "Edge from " << i << " to " << j << " is missing.";
      }
    }
  }
}

// Test case 4: neighbors - verify correct neighbor list for a given node.
TEST(PDAGTest, Neighbors) {
  constexpr size_t n = 5;
  PDAG g(n);

  // Initially, node 2 should have no neighbors.
  auto nbrs = g.neighbors(2);
  EXPECT_TRUE(nbrs.empty()) << "Node 2 should have no neighbors initially.";

  // Add edges from node 2 to nodes 0, 1, and 3.
  g.add_edge(2, 0);
  g.add_edge(2, 1);
  g.add_edge(2, 3);

  nbrs = g.neighbors(2);
  std::sort(nbrs.begin(), nbrs.end());
  std::vector<size_t> expected{0, 1, 3};
  EXPECT_EQ(nbrs, expected);

  // Node 4 still should have no neighbors.
  auto nbrs4 = g.neighbors(4);
  EXPECT_TRUE(nbrs4.empty());
}
