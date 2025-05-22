#include <digest/mod_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <digest/window_minimizer.hpp>
#include <digest/thread_out.hpp>
#include <variant>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

namespace py = pybind11;

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
window_minimizer_wrapper(const char *seq, size_t len, unsigned k, unsigned large_window,
				unsigned threads, bool include_hash = false) {
	if (threads == 1) {
		digest::WindowMin<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>
			digester(seq, len, k, large_window);
		if (include_hash) {
			std::vector<std::pair<uint32_t, uint32_t>> output;
			digester.roll_minimizer(len, output);
			return output;
		} else {
			std::vector<uint32_t> output;
			digester.roll_minimizer(len, output);
			return output;
		}
	}
	// return {};
	else {
		if (include_hash) {
			std::vector<std::vector<std::pair<uint32_t, uint32_t>>> output;
			digest::thread_out::thread_wind<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>(
				threads, output, seq, len, k, large_window);
			
			size_t total_size = 0;
			for (const auto& vec : output) {
				total_size += vec.size();
			}
			std::vector<std::pair<uint32_t, uint32_t>> concatenated;
			concatenated.reserve(total_size);
			for (const auto& vec : output) {
				concatenated.insert(concatenated.end(), vec.begin(), vec.end());
			}
			return concatenated;
		} else {
			std::vector<std::vector<uint32_t>> output;
			digest::thread_out::thread_wind<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>(
				threads, output, seq, len, k, large_window);
			
			size_t total_size = 0;
			for (const auto& vec : output) {
				total_size += vec.size();
			}
			
			std::vector<uint32_t> concatenated;
			concatenated.reserve(total_size);
		for (const auto& vec : output) {
				concatenated.insert(concatenated.end(), vec.begin(), vec.end());
			}
			return concatenated;
		}
	}
}

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
modimizer_wrapper(const char *seq, size_t len, unsigned k, uint32_t mod,
		  unsigned threads, bool include_hash = false) {
	if (threads == 1) {
		digest::ModMin<digest::BadCharPolicy::SKIPOVER> digester(seq, len, k, mod);
		if (include_hash) {
			std::vector<std::pair<uint32_t, uint32_t>> output;
			digester.roll_minimizer(len, output);
			return output;
		} else {
			std::vector<uint32_t> output;
			digester.roll_minimizer(len, output);
			return output;
		}
	}
	else {
		if (include_hash) {
			std::vector<std::vector<std::pair<uint32_t, uint32_t>>> output;
			digest::thread_out::thread_mod<digest::BadCharPolicy::SKIPOVER>(
				threads, output, seq, len, k, mod);
			
			size_t total_size = 0;
			for (const auto& vec : output) {
				total_size += vec.size();
			}
			std::vector<std::pair<uint32_t, uint32_t>> concatenated;
			concatenated.reserve(total_size);
			for (const auto& vec : output) {
				concatenated.insert(concatenated.end(), vec.begin(), vec.end());
			}
			return concatenated;
		} else {
			std::vector<std::vector<uint32_t>> output;
			digest::thread_out::thread_mod<digest::BadCharPolicy::SKIPOVER>(
				threads, output, seq, len, k, mod);
			
			size_t total_size = 0;
			for (const auto& vec : output) {
				total_size += vec.size();
			}
			
			std::vector<uint32_t> concatenated;
			concatenated.reserve(total_size);
			for (const auto& vec : output) {
				concatenated.insert(concatenated.end(), vec.begin(), vec.end());
			}
			return concatenated;
		}
	}
}

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
syncmer_wrapper(const char *seq, size_t len, unsigned k, unsigned large_window,
		unsigned threads, bool include_hash = false) {
	if (threads == 1) {
		digest::Syncmer<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>
			digester(seq, len, k, large_window);
		if (include_hash) {
			std::vector<std::pair<uint32_t, uint32_t>> output;
			digester.roll_minimizer(len, output);
			return output;
		} else {
			std::vector<uint32_t> output;
			digester.roll_minimizer(len, output);
			return output;
		}
	}
	else {
		if (include_hash) {
			std::vector<std::vector<std::pair<uint32_t, uint32_t>>> output;
			digest::thread_out::thread_sync<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>(
				threads, output, seq, len, k, large_window);
			
			size_t total_size = 0;
			for (const auto& vec : output) {
				total_size += vec.size();
			}
			std::vector<std::pair<uint32_t, uint32_t>> concatenated;
			concatenated.reserve(total_size);
			for (const auto& vec : output) {
				concatenated.insert(concatenated.end(), vec.begin(), vec.end());
			}
			return concatenated;
		} else {
			std::vector<std::vector<uint32_t>> output;
			digest::thread_out::thread_sync<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>(
				threads, output, seq, len, k, large_window);
			
			size_t total_size = 0;
			for (const auto& vec : output) {
				total_size += vec.size();
			}
			
			std::vector<uint32_t> concatenated;
			concatenated.reserve(total_size);
			for (const auto& vec : output) {
				concatenated.insert(concatenated.end(), vec.begin(), vec.end());
			}
			return concatenated;
		}
	}
}

pybind11::array window_minimizer_numpy(py::bytes seq, unsigned k, unsigned large_window,
                                       unsigned threads, bool include_hash = false) {
	py::gil_scoped_acquire gil;									
    py::buffer_info info(py::buffer(seq).request());
    const char* data = static_cast<const char*>(info.ptr);
    size_t len = static_cast<size_t>(info.size);

    std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>> result;
    {
        py::gil_scoped_release release;
        result = window_minimizer_wrapper(data, len, k, large_window, threads, include_hash);
    }
    if (std::holds_alternative<std::vector<uint32_t>>(result)) {
        const auto& vec = std::get<std::vector<uint32_t>>(result);
        return pybind11::array_t<uint32_t>(vec.size(), vec.data());
    } else {
        const auto& vec = std::get<std::vector<std::pair<uint32_t, uint32_t>>>(result);
        pybind11::array_t<uint32_t> arr({static_cast<pybind11::ssize_t>(vec.size()), static_cast<pybind11::ssize_t>(2)});
        auto buf = arr.mutable_unchecked<2>();
        for (pybind11::ssize_t i = 0; i < static_cast<pybind11::ssize_t>(vec.size()); ++i) {
            buf(i, 0) = vec[i].first;
            buf(i, 1) = vec[i].second;
        }
        return arr;
    }
}

pybind11::array modimizer_numpy(py::bytes seq, unsigned k, uint32_t mod,
                                unsigned threads, bool include_hash = false) {
	py::gil_scoped_acquire gil;
    py::buffer_info info(py::buffer(seq).request());
    const char* data = static_cast<const char*>(info.ptr);
    size_t len = static_cast<size_t>(info.size);

    std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>> result;
    {
        py::gil_scoped_release release;
        result = modimizer_wrapper(data, len, k, mod, threads, include_hash);
    }
    if (std::holds_alternative<std::vector<uint32_t>>(result)) {
        const auto& vec = std::get<std::vector<uint32_t>>(result);
        return pybind11::array_t<uint32_t>(vec.size(), vec.data());
    } else {
        const auto& vec = std::get<std::vector<std::pair<uint32_t, uint32_t>>>(result);
        pybind11::array_t<uint32_t> arr({static_cast<pybind11::ssize_t>(vec.size()), static_cast<pybind11::ssize_t>(2)});
        auto buf = arr.mutable_unchecked<2>();
        for (pybind11::ssize_t i = 0; i < static_cast<pybind11::ssize_t>(vec.size()); ++i) {
            buf(i, 0) = vec[i].first;
            buf(i, 1) = vec[i].second;
        }
        return arr;
    }
}

pybind11::array syncmer_numpy(py::bytes seq, unsigned k, unsigned large_window,
                              unsigned threads, bool include_hash = false) {
	py::gil_scoped_acquire gil;							
    py::buffer_info info(py::buffer(seq).request());
    const char* data = static_cast<const char*>(info.ptr);
    size_t len = static_cast<size_t>(info.size);

    std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>> result;
    {
        py::gil_scoped_release release;
        result = syncmer_wrapper(data, len, k, large_window, threads, include_hash);
    }
    if (std::holds_alternative<std::vector<uint32_t>>(result)) {
        const auto& vec = std::get<std::vector<uint32_t>>(result);
        return pybind11::array_t<uint32_t>(vec.size(), vec.data());
    } else {
        const auto& vec = std::get<std::vector<std::pair<uint32_t, uint32_t>>>(result);
        pybind11::array_t<uint32_t> arr({static_cast<pybind11::ssize_t>(vec.size()), static_cast<pybind11::ssize_t>(2)});
        auto buf = arr.mutable_unchecked<2>();
        for (pybind11::ssize_t i = 0; i < static_cast<pybind11::ssize_t>(vec.size()); ++i) {
            buf(i, 0) = vec[i].first;
            buf(i, 1) = vec[i].second;
        }
        return arr;
    }
}

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
window_minimizer(py::bytes seq, unsigned k, unsigned large_window,
                    unsigned threads, bool include_hash = false) {
    py::gil_scoped_acquire gil;
    py::buffer_info info(py::buffer(seq).request());
    const char* data = static_cast<const char*>(info.ptr);
    size_t len = static_cast<size_t>(info.size);
    {
        py::gil_scoped_release release;
        return window_minimizer_wrapper(data, len, k, large_window, threads, include_hash);
    }
}

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
modimizer(py::bytes seq, unsigned k, uint32_t mod,
             unsigned threads, bool include_hash = false) {
    py::gil_scoped_acquire gil;
    py::buffer_info info(py::buffer(seq).request());
    const char* data = static_cast<const char*>(info.ptr);
    size_t len = static_cast<size_t>(info.size);
    {
        py::gil_scoped_release release;
        return modimizer_wrapper(data, len, k, mod, threads, include_hash);
    }
}

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
syncmer(py::bytes seq, unsigned k, unsigned large_window,
           unsigned threads, bool include_hash = false) {
    py::gil_scoped_acquire gil;
    py::buffer_info info(py::buffer(seq).request());
    const char* data = static_cast<const char*>(info.ptr);
    size_t len = static_cast<size_t>(info.size);
    {
        py::gil_scoped_release release;
        return syncmer_wrapper(data, len, k, large_window, threads, include_hash);
    }
}
