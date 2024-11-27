#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

class CPDAGwithAdjMat {
private:
    size_t n;
    std::vector<uint64_t> adj_mat;

public:
    CPDAGwithAdjMat(size_t n);
    bool has_edge(int from, int to);
    add_edge(int from, int to);
    remove_edge(int from, int to);
    
    std::unordered_map<int, std::vector<int>> get_adj_list();
};