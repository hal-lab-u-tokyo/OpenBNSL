#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

/**
 * @class PDAGwithAdjMat
 * @brief Represents a Partially Directed Acyclic Graph (PDAG).
 *
 * This class implements a PDAG where edges are stored in a bit-packed adjacency matrix.
 * It is designed for Constraint-based Structure Learning algorithms.
 * This class supports exporting to `pgmpy` via an adjacency list representation.
 */
class PDAGwithAdjMat {
private:
    size_t n;                       ///< Number of nodes in the graph.
    std::vector<uint64_t> adj_mat;  ///< Bit-packed adjacency matrix.

public:
    /**
     * @brief Constructs a PDAGwithAdjMat object.
     * @param n The number of nodes in the graph.
     *
     * Initializes an adjacency matrix for a graph with `n` nodes.
     * All edges are initially set to 1 (i.e., complete undirected graph).
     */
    PDAGwithAdjMat(size_t n);

    /**
     * @brief Checks if there is a directed edge between two nodes.
     * @param from The source node.
     * @param to The destination node.
     * @return `true` if the edge exists, otherwise `false`.
     *
     * This method checks the adjacency matrix to determine if there
     * is a directed edge from `from` to `to`.
     */
    bool has_edge(int from, int to);

    /**
     * @brief Adds a directed edge between two nodes.
     * @param from The source node.
     * @param to The destination node.
     * @throws std::runtime_error If the edge already exists.
     *
     * This method sets the corresponding bit in the adjacency matrix
     * to indicate the presence of a directed edge from `from` to `to`.
     */
    void add_edge(int from, int to);

    /**
     * @brief Removes a directed edge between two nodes.
     * @param from The source node.
     * @param to The destination node.
     * @throws std::runtime_error If the edge does not exist.
     *
     * This method clears the corresponding bit in the adjacency matrix
     * to indicate the removal of a directed edge from `from` to `to`.
     */
    void remove_edge(int from, int to);
    
    /**
     * @brief Converts the graph into an adjacency list representation.
     * @return A map where each key is a node, and the value is a vector of nodes
     *         to which the key node has outgoing edges.
     *
     * This method traverses the adjacency matrix and constructs an
     * adjacency list, which is useful for integration with external
     * libraries such as `pgmpy`.
     */
    std::unordered_map<int, std::vector<int>> get_adj_list();
};