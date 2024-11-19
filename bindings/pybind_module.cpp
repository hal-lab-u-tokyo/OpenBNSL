#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "template.h"

namespace py = pybind11;

PYBIND11_MODULE(openbn, m) {
    m.doc() = "OpenBN Python bindings";

    // Bind the Template class
    py::class_<Template>(m, "Template")
        .def(py::init<const std::string&>())
        .def("print", &Template::print);

    
}
