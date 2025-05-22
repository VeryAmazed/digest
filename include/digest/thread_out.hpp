#ifndef THREAD_OUT_HPP
#define THREAD_OUT_HPP

#include "digest/mod_minimizer.hpp"
#include "digest/syncmer.hpp"
#include "digest/window_minimizer.hpp"
#include <cstdint>
#include <future>
#include <thread>
#include <vector>

/**
 *
 * @brief Possible implementation for multi-threading the digestion of a single
 * sequence. The key thing to note is basically by carefully telling where each
 * digester should start digesting you can make it so each kmer is only
 * considered once.
 * For more details on a function, click on more and it will take you to the
 * description that is located in modules
 *
 * @par IMPORTANT:
 * This approach will not generate correct results for sequences
 * that contain non-ACTG characters. Take this example, seq = ACTGANACNACTGA, k
 * = 4, l_wind = 4, thread_count = 2, there is a total of 4 valid kmers in this
 * sequence, and thus only 1 valid large window, but we can't know this until it
 * actually goes through the sequence, so it's going to try to partition the
 * sequence into ACTGANACNA, and ANACNACTGA and feed it into 2 digester objects
 * which now each have 0 valid large windows
 */
namespace digest::thread_out {

/**
 * @brief Exception thrown when invalid parameters are passed to the thread
 * functions
 */
class BadThreadOutParams : public std::exception {
	const char *what() const throw() {
		return "k must be greater than 3, start must be less than len, \
        and num threads must be greater or equal to the number of kmers/large windows \
        large_wind_kmer_am can't be 0";
	}
};

//------------- WORKER FUNCTIONS ----------------

// function that's passed to the thread for ModMinmizers
template <digest::BadCharPolicy P>
std::vector<uint32_t> thread_mod_roll1(const char *seq, size_t ind, unsigned k,
									   uint32_t mod, uint32_t congruence,
									   digest::MinimizedHashType minimized_h,
									   unsigned assigned_kmer_am) {
	std::vector<uint32_t> out;
	digest::ModMin<P> dig(seq, ind + assigned_kmer_am + k - 1, k, mod,
						  congruence, ind, minimized_h);
	dig.roll_minimizer(assigned_kmer_am, out);
	return out;
}

template <digest::BadCharPolicy P>
std::vector<std::pair<uint32_t, uint32_t>>
thread_mod_roll2(const char *seq, size_t ind, unsigned k, uint32_t mod,
				 uint32_t congruence, digest::MinimizedHashType minimized_h,
				 unsigned assigned_kmer_am) {
	std::vector<std::pair<uint32_t, uint32_t>> out;
	digest::ModMin<P> dig(seq, ind + assigned_kmer_am + k - 1, k, mod,
						  congruence, ind, minimized_h);
	dig.roll_minimizer(assigned_kmer_am, out);
	return out;
}

// function that's passed to the thread for WindowMinimizers
template <digest::BadCharPolicy P, class T>
std::vector<uint32_t> thread_wind_roll1(const char *seq, size_t ind, unsigned k,
										uint32_t large_wind_kmer_am,
										digest::MinimizedHashType minimized_h,
										unsigned assigned_lwind_am) {
	std::vector<uint32_t> out;
	digest::WindowMin<P, T> dig(
		seq, ind + assigned_lwind_am + k + large_wind_kmer_am - 1 - 1, k,
		large_wind_kmer_am, ind, minimized_h);
	dig.roll_minimizer(assigned_lwind_am, out);
	return out;
}

template <digest::BadCharPolicy P, class T>
std::vector<std::pair<uint32_t, uint32_t>> thread_wind_roll2(
	const char *seq, size_t ind, unsigned k, uint32_t large_wind_kmer_am,
	digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am) {
	std::vector<std::pair<uint32_t, uint32_t>> out;
	digest::WindowMin<P, T> dig(
		seq, ind + assigned_lwind_am + k + large_wind_kmer_am - 1 - 1, k,
		large_wind_kmer_am, ind, minimized_h);
	dig.roll_minimizer(assigned_lwind_am, out);
	return out;
}

// function that's passed to the thread for Syncmers
template <digest::BadCharPolicy P, class T>
std::vector<uint32_t> thread_sync_roll1(const char *seq, size_t ind, unsigned k,
										uint32_t large_wind_kmer_am,
										digest::MinimizedHashType minimized_h,
										unsigned assigned_lwind_am) {
	std::vector<uint32_t> out;
	digest::Syncmer<P, T> dig(
		seq, ind + assigned_lwind_am + k + large_wind_kmer_am - 1 - 1, k,
		large_wind_kmer_am, ind, minimized_h);
	dig.roll_minimizer(assigned_lwind_am, out);
	return out;
}

template <digest::BadCharPolicy P, class T>
std::vector<std::pair<uint32_t, uint32_t>> thread_sync_roll2(
	const char *seq, size_t ind, unsigned k, uint32_t large_wind_kmer_am,
	digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am) {
	std::vector<std::pair<uint32_t, uint32_t>> out;
	digest::Syncmer<P, T> dig(
		seq, ind + assigned_lwind_am + k + large_wind_kmer_am - 1 - 1, k,
		large_wind_kmer_am, ind, minimized_h);
	dig.roll_minimizer(assigned_lwind_am, out);
	return out;
}

/**
 * @param thread_count the number of threads to use
 * @param vec a vector of vectors in which the minimizers will be placed.
 *      Each vector corresponds to one thread. The minimizers within each vector
 *      will be in ascending order by index, and the vectors themselves will
 * also be in ascending order by index, i.e. all minimizers in vector_i will go
 *      before all minimizers in vector_(i+1).
 * @param seq char pointer poitning to the c-string of DNA sequence to be
 * hashed.
 * @param len length of seq.
 * @param k k-mer size.
 * @param mod mod space to be used to calculate universal minimizers
 * @param congruence value we want minimizer hashes to be congruent to in the
 * mod space
 * @param start 0-indexed position in seq to start hashing from.
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2
 * for reverse
 *
 * @throws BadThreadOutParams
 */
template <digest::BadCharPolicy P>
void thread_mod(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const char *seq, size_t len, unsigned k, uint32_t mod,
	uint32_t congruence = 0, size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	int num_kmers = (int)len - (int)start - (int)k + 1;
	if (k < 4 || start >= len || num_kmers < 0 ||
		(unsigned)num_kmers < thread_count) {
		throw BadThreadOutParams();
	}
	unsigned kmers_per_thread = num_kmers / thread_count;
	unsigned extras = num_kmers % thread_count;
	vec.reserve(thread_count);
	std::vector<std::future<std::vector<uint32_t>>> thread_vector;

	size_t ind = start;
	for (unsigned i = 0; i < thread_count; i++) {
		// issue is here
		// this will lead to a leak
		unsigned assigned_kmer_am = kmers_per_thread;
		if (extras > 0) {
			++(assigned_kmer_am);
			extras--;
		}

		thread_vector.emplace_back(std::async(thread_mod_roll1<P>, seq, ind, k,
											  mod, congruence, minimized_h,
											  assigned_kmer_am));

		ind += assigned_kmer_am;
	}
	for (auto &t : thread_vector) {
		vec.emplace_back(std::move(t.get()));
	}
}

/**
 * @brief same as the other thread_mod, except it can take a C++ string, and
 * does not need to be provided the length of the string
 *
 * @param seq C++ string of DNA sequence to be hashed.
 */
template <digest::BadCharPolicy P>
void thread_mod(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const std::string &seq, unsigned k, uint32_t mod, uint32_t congruence = 0,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	thread_mod<P>(thread_count, vec, seq.c_str(), seq.size(), k, mod,
				  congruence, start, minimized_h);
}

/**
 * @brief same as other thread_mod that takes a c-string,
 * except here vec is a vector of vectors of pairs of uint32_ts
 *
 * @param vec vec will contain both the index and the hash of minimizers.
 * All other things previously stated about vec remain true
 */
template <digest::BadCharPolicy P>
void thread_mod(
	unsigned thread_count,
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &vec,
	const char *seq, size_t len, unsigned k, uint32_t mod,
	uint32_t congruence = 0, size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	int num_kmers = (int)len - (int)start - (int)k + 1;
	if (k < 4 || start >= len || num_kmers < 0 ||
		(unsigned)num_kmers < thread_count) {
		throw BadThreadOutParams();
	}
	unsigned kmers_per_thread = num_kmers / thread_count;
	unsigned extras = num_kmers % thread_count;
	vec.reserve(thread_count);
	std::vector<std::future<std::vector<std::pair<uint32_t, uint32_t>>>>
		thread_vector;

	size_t ind = start;
	for (unsigned i = 0; i < thread_count; i++) {
		// issue is here
		// this will lead to a leak
		unsigned assigned_kmer_am = kmers_per_thread;
		if (extras > 0) {
			++(assigned_kmer_am);
			extras--;
		}

		thread_vector.emplace_back(std::async(thread_mod_roll2<P>, seq, ind, k,
											  mod, congruence, minimized_h,
											  assigned_kmer_am));

		ind += assigned_kmer_am;
	}
	for (auto &t : thread_vector) {
		vec.emplace_back(std::move(t.get()));
	}
}

/**
 * @brief same as other thread_mod that takes a C++ string,
 * except here vec is a vector of vectors of pairs of uint32_ts
 *
 * @param vec vec will contain both the index and the hash of minimizers.
 * All other things previously stated about vec remain true
 */
template <digest::BadCharPolicy P>
void thread_mod(
	unsigned thread_count,
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &vec,
	const std::string &seq, unsigned k, uint32_t mod, uint32_t congruence = 0,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	thread_mod<P>(thread_count, vec, seq.c_str(), seq.size(), k, mod,
				  congruence, start, minimized_h);
}

/**
 * @tparam P policy for dealing with non-ACTG characters
 * @tparam T min query data structure to use, refer to docs of the classes in
 * the ds namespace for more info
 *
 * @param thread_count the number of threads to use
 * @param vec a vector of vectors in which the minimizers will be placed.
 *      Each vector corresponds to one thread. The minimizers within each vector
 *      will be in ascending order by index, and the vectors themselves will
 * also be in ascending order by index, i.e. all minimizers in vector_i will go
 *      before all minimizers in vector_(i+1).
 * @param seq char pointer poitning to the c-string of DNA sequence to be
 * hashed.
 * @param len length of seq.
 * @param k k-mer size.
 * @param large_wind_kmer_am
 * @param start 0-indexed position in seq to start hashing from.
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2
 * for reverse
 *
 * @throws BadThreadOutParams
 */
template <digest::BadCharPolicy P, class T>
void thread_wind(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const char *seq, size_t len, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	int num_lwinds = (int)len - (int)start - (int)(k + large_wind_kmer_am) + 2;
	if (large_wind_kmer_am == 0 || k < 4 || start >= len || num_lwinds < 0 ||
		(unsigned)num_lwinds < thread_count) {
		throw BadThreadOutParams();
	}
	unsigned lwinds_per_thread = num_lwinds / thread_count;
	unsigned extras = num_lwinds % thread_count;
	vec.reserve(thread_count);
	std::vector<std::future<std::vector<uint32_t>>> thread_vector;

	size_t ind = start;
	for (unsigned i = 0; i < thread_count; i++) {
		// issue is here
		// this will lead to a leak
		unsigned assigned_lwind_am = lwinds_per_thread;
		if (extras > 0) {
			++(assigned_lwind_am);
			extras--;
		}

		thread_vector.push_back(std::async(thread_wind_roll1<P, T>, seq, ind,
											  k, large_wind_kmer_am,
											  minimized_h, assigned_lwind_am));

		ind += assigned_lwind_am;
	}
	for (auto &t : thread_vector) {
		vec.emplace_back(std::move(t.get()));
	}

	// handle duplicates
	// the only possible place for a duplicate is for the last element
	// of vec[i] to equal the first value of vec[i+1] due to the fact
	// that thread_i+1 can't know the last minimizer of thread_i
	for (unsigned i = 0; i < thread_count - 1; i++) {
		int last = (int)vec[i].size() - 1;
		if (vec[i][last] == vec[i + 1][0]) {
			vec[i].pop_back();
		}
	}
}

/**
 * @brief same as the other thread_wind, except it can take a C++ string, and
 * does not need to be provided the length of the string
 *
 * @param seq C++ string of DNA sequence to be hashed.
 */
template <digest::BadCharPolicy P, class T>
void thread_wind(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const std::string &seq, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	thread_wind<P, T>(thread_count, vec, seq.c_str(), seq.size(), k,
					  large_wind_kmer_am, start, minimized_h);
}

/**
 * @brief same as other thread_wind that takes a c-string,
 * except here vec is a vector of vectors of pairs of uint32_ts
 *
 * @param vec vec will contain both the index and the hash of minimizers.
 * All other things previously stated about vec remain true
 */
template <digest::BadCharPolicy P, class T>
void thread_wind(
	unsigned thread_count,
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &vec,
	const char *seq, size_t len, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	int num_lwinds = (int)len - (int)start - (int)(k + large_wind_kmer_am) + 2;
	if (large_wind_kmer_am == 0 || k < 4 || start >= len || num_lwinds < 0 ||
		(unsigned)num_lwinds < thread_count) {
		throw BadThreadOutParams();
	}
	unsigned lwinds_per_thread = num_lwinds / thread_count;
	unsigned extras = num_lwinds % thread_count;
	vec.reserve(thread_count);
	std::vector<std::future<std::vector<std::pair<uint32_t, uint32_t>>>> thread_vector;

	size_t ind = start;
	for (unsigned i = 0; i < thread_count; i++) {
		// issue is here
		// this will lead to a leak
		unsigned assigned_lwind_am = lwinds_per_thread;
		if (extras > 0) {
			++(assigned_lwind_am);
			extras--;
		}

		thread_vector.push_back(std::async(thread_wind_roll2<P, T>, seq, ind,
											  k, large_wind_kmer_am,
											  minimized_h, assigned_lwind_am));

		ind += assigned_lwind_am;
	}
	for (auto &t : thread_vector) {
		vec.emplace_back(std::move(t.get()));
	}

	// handle duplicates
	// the only possible place for a duplicate is for the last element
	// of vec[i] to equal the first value of vec[i+1] due to the fact
	// that thread_i+1 can't know the last minimizer of thread_i
	for (unsigned i = 0; i < thread_count - 1; i++) {
		int last = (int)vec[i].size() - 1;
		if (vec[i][last] == vec[i + 1][0]) {
			vec[i].pop_back();
		}
	}
}

/**
 * @brief same as other thread_wind that takes a C++ string,
 * except here vec is a vector of vectors of pairs of uint32_ts
 *
 * @param vec vec will contain both the index and the hash of minimizers.
 * All other things previously stated about vec remain true
 */
template <digest::BadCharPolicy P, class T>
void thread_wind(
	unsigned thread_count,
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &vec,
	const std::string &seq, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	thread_wind<P, T>(thread_count, vec, seq.c_str(), seq.size(), k,
					  large_wind_kmer_am, start, minimized_h);
}

/**
 * @tparam P policy for dealing with non-ACTG characters
 * @tparam T min query data structure to use, refer to docs of the classes in
 * the ds namespace for more info
 *
 * @param thread_count the number of threads to use
 * @param vec a vector of vectors in which the minimizers will be placed.
 *      Each vector corresponds to one thread. The minimizers within each vector
 *      will be in ascending order by index, and the vectors themselves will
 * also be in ascending order by index, i.e. all minimizers in vector_i will go
 *      before all minimizers in vector_(i+1).
 * @param seq char pointer poitning to the c-string of DNA sequence to be
 * hashed.
 * @param len length of seq.
 * @param k k-mer size.
 * @param large_wind_kmer_am
 * @param start 0-indexed position in seq to start hashing from.
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2
 * for reverse
 *
 * @throws BadThreadOutParams
 */
template <digest::BadCharPolicy P, class T>
void thread_sync(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const char *seq, size_t len, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	int num_lwinds = (int)len - (int)start - (int)(k + large_wind_kmer_am) + 2;
	if (large_wind_kmer_am == 0 || k < 4 || start >= len || num_lwinds < 0 ||
		(unsigned)num_lwinds < thread_count) {
		throw BadThreadOutParams();
	}
	unsigned lwinds_per_thread = num_lwinds / thread_count;
	unsigned extras = num_lwinds % thread_count;
	vec.reserve(thread_count);
	std::vector<std::future<std::vector<uint32_t>>> thread_vector;

	size_t ind = start;
	for (unsigned i = 0; i < thread_count; i++) {
		// issue is here
		// this will lead to a leak
		unsigned assigned_lwind_am = lwinds_per_thread;
		if (extras > 0) {
			++(assigned_lwind_am);
			extras--;
		}

		thread_vector.emplace_back(std::async(thread_sync_roll1<P, T>, seq, ind,
											  k, large_wind_kmer_am,
											  minimized_h, assigned_lwind_am));

		ind += assigned_lwind_am;
	}
	for (auto &t : thread_vector) {

		vec.emplace_back(std::move(t.get()));
	}
}

/**
 * @brief same as the other thread_sync, except it can take a C++ string, and
 * does not need to be provided the length of the string
 *
 * @param seq C++ string of DNA sequence to be hashed.
 */
template <digest::BadCharPolicy P, class T>
void thread_sync(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const std::string &seq, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	thread_sync<P, T>(thread_count, vec, seq.c_str(), seq.size(), k,
					  large_wind_kmer_am, start, minimized_h);
}

/**
 * @brief same as other thread_wind that takes a c-string,
 * except here vec is a vector of vectors of pairs of uint32_ts
 *
 * @param vec vec will contain both the index and the hash of minimizers.
 * All other things previously stated about vec remain true
 */
template <digest::BadCharPolicy P, class T>
void thread_sync(
	unsigned thread_count,
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &vec,
	const char *seq, size_t len, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	int num_lwinds = (int)len - (int)start - (int)(k + large_wind_kmer_am) + 2;
	if (large_wind_kmer_am == 0 || k < 4 || start >= len || num_lwinds < 0 ||
		(unsigned)num_lwinds < thread_count) {
		throw BadThreadOutParams();
	}
	unsigned lwinds_per_thread = num_lwinds / thread_count;
	unsigned extras = num_lwinds % thread_count;
	vec.reserve(thread_count);
	std::vector<std::future<std::vector<std::pair<uint32_t, uint32_t>>>>
		thread_vector;

	size_t ind = start;
	for (unsigned i = 0; i < thread_count; i++) {
		// issue is here
		// this will lead to a leak
		unsigned assigned_lwind_am = lwinds_per_thread;
		if (extras > 0) {
			++(assigned_lwind_am);
			extras--;
		}

		thread_vector.emplace_back(std::async(thread_sync_roll2<P, T>, seq, ind,
											  k, large_wind_kmer_am,
											  minimized_h, assigned_lwind_am));

		ind += assigned_lwind_am;
	}
	for (auto &t : thread_vector) {
		vec.emplace_back(std::move(t.get()));
	}
}

/**
 * @brief same as other thread_sync that takes a C++ string,
 * except here vec is a vector of vectors of pairs of uint32_ts
 *
 * @param vec vec will contain both the index and the hash of minimizers.
 * All other things previously stated about vec remain true
 */
template <digest::BadCharPolicy P, class T>
void thread_sync(
	unsigned thread_count,
	std::vector<std::vector<std::pair<uint32_t, uint32_t>>> &vec,
	const std::string &seq, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON) {
	thread_sync<P, T>(thread_count, vec, seq.c_str(), seq.size(), k,
					  large_wind_kmer_am, start, minimized_h);
}

} // namespace digest::thread_out

#endif // THREAD_OUT_HPP
