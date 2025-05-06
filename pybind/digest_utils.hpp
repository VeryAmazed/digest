#include <digest/mod_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <digest/window_minimizer.hpp>
#include <digest/thread_out.hpp>
#include <variant>

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>>
window_minimizer(const std::string &seq, unsigned k, unsigned large_window,
				unsigned threads, bool include_hash = false) {
	if (threads == 1) {
		digest::WindowMin<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>
			digester(seq, k, large_window);
		if (include_hash) {
			std::vector<std::pair<uint32_t, uint32_t>> output;
			digester.roll_minimizer(seq.length(), output);
			return output;
		} else {
			std::vector<uint32_t> output;
			digester.roll_minimizer(seq.length(), output);
			return output;
		}
	}
	// return {};
	else {
		if (include_hash) {
			std::vector<std::vector<std::pair<uint32_t, uint32_t>>> output;
			digest::thread_out::thread_wind<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>(
				threads, output, seq, k, large_window);
			
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
				threads, output, seq, k, large_window);
			
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
modimizer(const std::string &seq, unsigned k, uint32_t mod,
		  unsigned threads, bool include_hash = false) {
	if (threads == 1) {
		digest::ModMin<digest::BadCharPolicy::SKIPOVER> digester(seq, k, mod);
		if (include_hash) {
			std::vector<std::pair<uint32_t, uint32_t>> output;
			digester.roll_minimizer(seq.length(), output);
			return output;
		} else {
			std::vector<uint32_t> output;
			digester.roll_minimizer(seq.length(), output);
			return output;
		}
	}
	else {
		if (include_hash) {
			std::vector<std::vector<std::pair<uint32_t, uint32_t>>> output;
			digest::thread_out::thread_mod<digest::BadCharPolicy::SKIPOVER>(
				threads, output, seq, k, mod);
			
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
				threads, output, seq, k, mod);
			
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
syncmer(const std::string &seq, unsigned k, unsigned large_window,
		unsigned threads, bool include_hash = false) {
	if (threads == 1) {
		digest::Syncmer<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>
			digester(seq, k, large_window);
		if (include_hash) {
			std::vector<std::pair<uint32_t, uint32_t>> output;
			digester.roll_minimizer(seq.length(), output);
			return output;
		} else {
			std::vector<uint32_t> output;
			digester.roll_minimizer(seq.length(), output);
			return output;
		}
	}
	else {
		if (include_hash) {
			std::vector<std::vector<std::pair<uint32_t, uint32_t>>> output;
			digest::thread_out::thread_sync<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive>(
				threads, output, seq, k, large_window);
			
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
				threads, output, seq, k, large_window);
			
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
