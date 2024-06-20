#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <digest_utils.hpp>

namespace py = pybind11;

PYBIND11_MODULE(Digest, m) {
    m.doc() = "bindings for Digest";
    m.def("window_minimizer", &window_minimizer, "A function that runs window minimizer digestion",
            py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11);
    m.def("modimizer", &modimizer, "A function that runs mod-minimizer digestion",
            py::arg("seq"), py::arg("k") = 31, py::arg("mod") = 100);
    m.def("syncmer", &syncmer, "A function that runs syncmer digestion",
            py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11);
}
