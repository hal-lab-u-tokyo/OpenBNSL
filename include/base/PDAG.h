#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

class PDAGwithAdjMat {
private:
    size_t n;
    std::vector<uint64_t> adj_mat;

public:
    PDAGwithAdjMat(size_t n);
    bool has_edge(int from, int to);
    void add_edge(int from, int to);
    void remove_edge(int from, int to);
    
    std::unordered_map<int, std::vector<int>> get_adj_list();
};