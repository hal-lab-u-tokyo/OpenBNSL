#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace py = pybind11;
#define N (10)

void bind_type_inspection(py::module& m) {
  auto submodule =
      m.def_submodule("type_inspection", "Type inspection submodule");

  m.def("test_bool", [](bool x) { return x; });
  m.def("test_int", [](int x) { return x; });
  m.def("test_unsigned_int", [](unsigned int x) { return x; });
  m.def("test_long", [](long x) { return x; });
  m.def("test_unsigned_long", [](unsigned long x) { return x; });
  m.def("test_float", [](float x) { return x; });
  m.def("test_double", [](double x) { return x; });
  m.def("test_string", [](const std::string& x) { return x; });
  m.def("test_pair_int", [](const std::pair<int, int>& x) { return x; });
  m.def("test_tuple_int", [](const std::tuple<int, int, int>& x) { return x; });
  m.def("test_vector_int", [](const std::vector<int>& x) { return x; });
  m.def("test_array_int", [](const std::array<int, N>& x) { return x; });
  m.def("test_map_string_int",
        [](const std::map<std::string, int>& x) { return x; });
  m.def("test_unordered_map_string_int",
        [](const std::unordered_map<std::string, int>& x) { return x; });
  m.def("test_set_int", [](const std::set<int>& x) { return x; });
  m.def("test_function_int_int",
        [](const std::function<int(int)>& x) { return x; });
}
