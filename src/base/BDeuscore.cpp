
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>
#include "StructureScore.cpp"  // Assuming BDeuScore is defined here

namespace py = pybind11;
using namespace Eigen;

MatrixXd pandas_to_eigen(const py::object& df) {
    auto values = df.attr("values").cast<py::array>();
    py::buffer_info info = values.request();
    double* ptr = static_cast<double*>(info.ptr);
    return Map<MatrixXd>(ptr, info.shape[0], info.shape[1]);
}

double BDeu_score(const py::object& df, int equivalent_sample_size) {
    MatrixXd data = pandas_to_eigen(df);
    return data;
    // BDeuScore BDeu_score(data, equivalent_sample_size);
    // BayesianNetwork model;  // Assuming BayesianNetwork is defined and can be constructed
    // // ...construct the BayesianNetwork model...
    // return BDeu_score.score(model);
}
