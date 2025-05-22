#include <digest_utils.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(digest, m) {
	m.doc() = "bindings for digest";
	
	m.def("window_minimizer", 
		  [](std::string_view seq, unsigned k, unsigned w, unsigned num_threads, bool include_hash) {
			  return window_minimizer_wrapper(seq.data(), seq.size(), k, w, num_threads, include_hash);
		  },
		  "A function that runs window minimizer digestion",
		  py::call_guard<py::gil_scoped_release>(),
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);

	m.def("modimizer", 
		  [](std::string_view seq, unsigned k, uint32_t mod, unsigned num_threads, bool include_hash) {
			  return modimizer_wrapper(seq.data(), seq.size(), k, mod, num_threads, include_hash);
		  },
		  "A function that runs mod-minimizer digestion",
		  py::call_guard<py::gil_scoped_release>(),
		  py::arg("seq"), py::arg("k") = 31, py::arg("mod") = 100, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);

	m.def("syncmer", 
		  [](std::string_view seq, unsigned k, unsigned w, unsigned num_threads, bool include_hash) {
			  return syncmer_wrapper(seq.data(), seq.size(), k, w, num_threads, include_hash);
		  },
		  "A function that runs syncmer digestion",
		  py::call_guard<py::gil_scoped_release>(),
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);

	m.def("window_minimizer_np",
		  [](std::string_view seq, unsigned k, unsigned w, unsigned num_threads, bool include_hash) {
			  return window_minimizer_numpy(seq.data(), seq.size(), k, w, num_threads, include_hash);
		  },
		  "A function that runs window minimizer digestion and returns a numpy array",
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("modimizer_np",
		  [](std::string_view seq, unsigned k, uint32_t mod, unsigned num_threads, bool include_hash) {
			  return modimizer_numpy(seq.data(), seq.size(), k, mod, num_threads, include_hash);
		  },
		  "A function that runs mod-minimizer digestion and returns a numpy array",
		  py::arg("seq"), py::arg("k") = 31, py::arg("mod") = 100, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
	m.def("syncmer_np",
		  [](std::string_view seq, unsigned k, unsigned w, unsigned num_threads, bool include_hash) {
			  return syncmer_numpy(seq.data(), seq.size(), k, w, num_threads, include_hash);
		  },
		  "A function that runs syncmer digestion and returns a numpy array",
		  py::arg("seq"), py::arg("k") = 31, py::arg("w") = 11, py::arg("num_threads") = 1,
		  py::arg("include_hash") = false);
}
