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

// Test case 5: CompleteGraph with node count as a multiple of 64 (e.g., n = 64)
TEST(PDAGTest, CompleteGraph_MultipleOf64) {
  constexpr size_t n = 64;  // num_vars is exactly a multiple of 64.
  PDAG g(n);
  g.complete_graph();

  // Self-loops must be removed.
  for (size_t i = 0; i < n; ++i) {
    EXPECT_FALSE(g.has_edge(i, i))
        << "Self-loop at node " << i << " should be removed.";
  }
  // Every edge from i to j (i != j) should be present.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i != j) {
        EXPECT_TRUE(g.has_edge(i, j))
            << "Edge from " << i << " to " << j << " is missing.";
      }
    }
  }
}

// Test case 6: CompleteGraph with node count not a multiple of 64 (e.g., n =
// 65)
TEST(PDAGTest, CompleteGraph_NonMultipleOf64) {
  constexpr size_t n = 65;  // num_vars is not a multiple of 64.
  PDAG g(n);
  g.complete_graph();

  // Self-loops must be removed.
  for (size_t i = 0; i < n; ++i) {
    EXPECT_FALSE(g.has_edge(i, i))
        << "Self-loop at node " << i << " should be removed.";
  }
  // Every edge from i to j (i != j) should be present.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i != j) {
        EXPECT_TRUE(g.has_edge(i, j))
            << "Edge from " << i << " to " << j << " is missing.";
      }
    }
  }
}

// Test case 7: Copy Constructor - Verify that the copy contains the same edges.
TEST(PDAGTest, CopyConstructor) {
  constexpr size_t n = 10;
  PDAG g(n);
  // Setup the original graph
  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(3, 4);
  g.add_edge(5, 9);

  PDAG g_copy(g);

  // Ensure all edges are identically copied.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      EXPECT_EQ(g.has_edge(i, j), g_copy.has_edge(i, j))
          << "Edge mismatch between original and copied graph at (" << i << ", "
          << j << ")";
    }
  }

  // Changing the original graph should not affect the copy.
  g.add_edge(2, 7);
  EXPECT_FALSE(g_copy.has_edge(2, 7))
      << "Copied graph should not change when the original is modified.";
}

// Test case 8: Assignment Operator - Verify that the assigned object has the
// same state.
TEST(PDAGTest, AssignmentOperator) {
  constexpr size_t n = 10;
  PDAG g1(n);
  PDAG g2(n);

  // Setup g1 with some edges.
  g1.add_edge(1, 2);
  g1.add_edge(2, 3);

  // Setup g2 with different edges (to be overwritten).
  g2.add_edge(5, 6);

  // Perform assignment.
  g2 = g1;

  // Verify that g2 now matches g1.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      EXPECT_EQ(g1.has_edge(i, j), g2.has_edge(i, j))
          << "After assignment, edge mismatch at (" << i << ", " << j << ")";
    }
  }

  // Further modifications to g1 should not affect g2.
  g1.add_edge(3, 4);
  EXPECT_FALSE(g2.has_edge(3, 4))
      << "Modification to original should not reflect in the assigned object.";
}

// Test case 9: Neighbors - Out-of-range index should throw an exception.
TEST(PDAGTest, Neighbors_OutOfRange) {
  constexpr size_t n = 5;
  PDAG g(n);

  EXPECT_THROW(g.neighbors(n), std::out_of_range);
}

// Test case 10: Directed and Undirected Edge Queries.
TEST(PDAGTest, DirectedUndirectedEdges) {
  constexpr size_t n = 4;
  PDAG g(n);

  // Add edge 0 -> 1 (directed).
  g.add_edge(0, 1);
  EXPECT_TRUE(g.has_edge(0, 1));
  EXPECT_FALSE(g.has_edge(1, 0));
  EXPECT_TRUE(g.has_directed_edge(0, 1));
  EXPECT_FALSE(g.has_directed_edge(1, 0));
  EXPECT_FALSE(g.has_undirected_edge(0, 1));

  // Add edge 1 -> 0 to form an undirected edge.
  g.add_edge(1, 0);
  EXPECT_TRUE(g.has_edge(1, 0));
  // Now, both directions exist; directed edge query should fail.
  EXPECT_FALSE(g.has_directed_edge(0, 1));
  // undirected edge holds.
  EXPECT_TRUE(g.has_undirected_edge(0, 1));
}

// Test case 11: Successors and Predecessors Queries.
TEST(PDAGTest, SuccessorsPredecessors) {
  constexpr size_t n = 6;
  PDAG g(n);

  // Construct the following graph:
  // 0 -> 1, 0 -> 2, 2 -> 3, 4 -> 3, 5 -> 3
  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(2, 3);
  g.add_edge(4, 3);
  g.add_edge(5, 3);

  // successors(0) should be {1, 2}.
  auto succ0 = g.successors(0);
  std::sort(succ0.begin(), succ0.end());
  std::vector<size_t> expected_succ0 = {1, 2};
  EXPECT_EQ(succ0, expected_succ0);

  // predecessors(3) should be {2, 4, 5}.
  auto pred3 = g.predecessors(3);
  std::sort(pred3.begin(), pred3.end());
  std::vector<size_t> expected_pred3 = {2, 4, 5};
  EXPECT_EQ(pred3, expected_pred3);
}

// Test case 12: Undirected Neighbors Query.
TEST(PDAGTest, UndirectedNeighbors) {
  constexpr size_t n = 4;
  PDAG g(n);

  // Create edges such that:
  // 0 <-> 1, 0 -> 2, and 3 -> 0.
  g.add_edge(0, 1);
  g.add_edge(1, 0);  // now 0 and 1 are connected undirectedly.
  g.add_edge(0, 2);  // only 0 -> 2.
  g.add_edge(3, 0);  // only 3 -> 0.

  // undirected_neighbors(0) should return only {1}.
  auto undir0 = g.undirected_neighbors(0);
  std::sort(undir0.begin(), undir0.end());
  std::vector<size_t> expected_undir0 = {1};
  EXPECT_EQ(undir0, expected_undir0);

  // For node 1, undirected_neighbors should return {0}.
  auto undir1 = g.undirected_neighbors(1);
  std::sort(undir1.begin(), undir1.end());
  std::vector<size_t> expected_undir1 = {0};
  EXPECT_EQ(undir1, expected_undir1);
}