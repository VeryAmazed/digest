#include <digest_utils.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(digest, m) {
	m.doc() = "bindings for digest";
	m.def("window_minimizer", &window_minimizer,
		  "A function that runs window minimizer digestion", 
		  py::call_guard<py::gil_scoped_release>(),  // Release GIL during execution
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("modimizer", &modimizer,
		  "A function that runs mod-minimizer digestion", 
		  py::call_guard<py::gil_scoped_release>(),  // Release GIL during execution
		  py::arg("seq"), py::arg("k") = 31, py::arg("mod") = 100, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("syncmer", &syncmer, "A function that runs syncmer digestion",
			py::call_guard<py::gil_scoped_release>(),  // Release GIL during execution
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("window_minimizer_numpy", &window_minimizer_numpy,
		  "A function that runs window minimizer digestion and returns a numpy array",
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("modimizer_numpy", &modimizer_numpy,
		  "A function that runs mod-minimizer digestion and returns a numpy array",
		  py::arg("seq"), py::arg("k") = 31, py::arg("mod") = 100, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("syncmer_numpy", &syncmer_numpy,
		  "A function that runs syncmer digestion and returns a numpy array",
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
}
