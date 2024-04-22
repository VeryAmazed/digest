#ifndef THREAD_OUT_HPP
#define THREAD_OUT_HPP

#include "digest/mod_minimizer.hpp"
#include "digest/syncmer.hpp"
#include "digest/window_minimizer.hpp"
#include <cstdint>
#include <thread>
#include <vector>

/**
 * \defgroup thread_out Threading functions
 * @{
 */

/**
 * \ingroup thread_out
 * 
 * @brief Description:
 * Possible implementation for multi-threading the digestion of a single
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

namespace thread_out {

class BadThreadOutParams : public std::exception {
	const char *what() const throw() {
		return "k must be greater than 3, start must be less than len, \
        and num threads must be greater or equal to the number of kmers/large windows \
        large_wind_kmer_am can't be 0";
	}
};

/**
 * \ingroup thread_out
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
 * @brief same as the other thread_mod, except it can take a C++ string, and does not need to be provided
 * the length of the string
 * 
 * @param seq C++ string of DNA sequence to be hashed.
 */
template <digest::BadCharPolicy P>
void thread_mod(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const std::string &seq, unsigned k, uint32_t mod, uint32_t congruence = 0,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
 * @tparam P policy for dealing with non-ACTG characters
 * @tparam T min query data structure to use, refer to docs of the classes in the ds namespace for more info
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
 * @brief same as the other thread_wind, except it can take a C++ string, and does not need to be provided
 * the length of the string
 * 
 * @param seq C++ string of DNA sequence to be hashed.
 */
template <digest::BadCharPolicy P, class T>
void thread_wind(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const std::string &seq, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 *  
 * @tparam P policy for dealing with non-ACTG characters
 * @tparam T min query data structure to use, refer to docs of the classes in the ds namespace for more info
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
 * @brief same as the other thread_sync, except it can take a C++ string, and does not need to be provided
 * the length of the string
 * 
 * @param seq C++ string of DNA sequence to be hashed.
 */
template <digest::BadCharPolicy P, class T>
void thread_sync(
	unsigned thread_count, std::vector<std::vector<uint32_t>> &vec,
	const std::string &seq, unsigned k, uint32_t large_wind_kmer_am,
	size_t start = 0,
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * \ingroup thread_out
 * 
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
	digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

//------------- WORKER FUNCTIONS ----------------

// function that's passed to the thread for ModMinmizers
template <digest::BadCharPolicy P>
std::vector<uint32_t> thread_mod_roll1(const char *seq, size_t ind, unsigned k,
									   uint32_t mod, uint32_t congruence,
									   digest::MinimizedHashType minimized_h,
									   unsigned assigned_kmer_am);

template <digest::BadCharPolicy P>
std::vector<std::pair<uint32_t, uint32_t>>
thread_mod_roll2(const char *seq, size_t ind, unsigned k, uint32_t mod,
				 uint32_t congruence, digest::MinimizedHashType minimized_h,
				 unsigned assigned_kmer_am);

// function that's passed to the thread for WindowMinimizers
template <digest::BadCharPolicy P, class T>
std::vector<uint32_t> thread_wind_roll1(const char *seq, size_t ind, unsigned k,
										uint32_t large_wind_kmer_am,
										digest::MinimizedHashType minimized_h,
										unsigned assigned_lwind_am);

template <digest::BadCharPolicy P, class T>
std::vector<std::pair<uint32_t, uint32_t>> thread_wind_roll2(
	const char *seq, size_t ind, unsigned k, uint32_t large_wind_kmer_am,
	digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am);

// function that's passed to the thread for Syncmers
template <digest::BadCharPolicy P, class T>
std::vector<uint32_t> thread_sync_roll1(const char *seq, size_t ind, unsigned k,
										uint32_t large_wind_kmer_am,
										digest::MinimizedHashType minimized_h,
										unsigned assigned_lwind_am);

template <digest::BadCharPolicy P, class T>
std::vector<std::pair<uint32_t, uint32_t>> thread_sync_roll2(
	const char *seq, size_t ind, unsigned k, uint32_t large_wind_kmer_am,
	digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am);

} // namespace thread_out
/**@}*/
#include "thread_out.tpp"
#endif
