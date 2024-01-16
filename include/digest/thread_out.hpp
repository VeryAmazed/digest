#ifndef THREAD_OUT_HPP
#define THREAD_OUT_HPP

#include "digest/mod_minimizer.hpp"
#include "digest/window_minimizer.hpp"
#include "digest/syncmer.hpp"
#include <cstdint>
#include <thread>

/*
    Possible implementation for multi-threading the digestion of a single sequence.
    The key thing to note is basically by carefully telling where each digester should start digesting
    you can make it so each kmer is only considered once.
    I have very little experience with threading, so you could probably thread this out better than me

    IMPORTANT: This approach will not generate correct results for sequences that contain non-ACTG
    characters. 
    Take this example, seq = ACTGANACNACTGA, k = 4, l_wind = 4, thread_count = 2, there
    is a total of 4 valid kmers in this sequence, and thus only 1 valid large window, but we 
    can't know this until it actually goes through the sequence, so it's going to try to partition the
    sequence into ACTGANACNA, and ANACNACTGA and feed it into 2 digester objects which now each have 0
    valid large windows
*/
namespace thread_out{

class BadThreadOutParams : public std::exception
{
	const char * what () const throw ()
    {
    	return "k must be greater than 3, start must be less than len, \
        and num threads must be greater or equal to the number of kmers/large windows \
        large_wind_kmer_am can't be 0";
    }
};

/**
 * @param thread_count the number of threads to use
 * @param vec a vector of vectors in which the minimizers will be placed.
 *      Each vector corresponds to one thread. The minimizers within each vector
 *      will be in ascending order by index, and the vectors themselves will also
 *      be in ascending order by index, i.e. all minimizers in vector_i will go
 *      before all minimizers in vector_(i+1). 
 * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
 * @param len length of seq.
 * @param k k-mer size.
 * @param mod mod space to be used to calculate universal minimizers
 * @param congruence value we want minimizer hashes to be congruent to in the mod space
 * @param start 0-indexed position in seq to start hashing from. 
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
 * 
 * @throws BadThreadOutParams 
 */
void thread_mod(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t mod, uint32_t congruence = 0, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
 */
void thread_mod(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t mod, uint32_t congruence = 0, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * @param thread_count the number of threads to use
 * @param vec a vector of vectors in which the minimizers will be placed.
 *      Each vector corresponds to one thread. The minimizers within each vector
 *      will be in ascending order by index, and the vectors themselves will also
 *      be in ascending order by index, i.e. all minimizers in vector_i will go
 *      before all minimizers in vector_(i+1). 
 * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
 * @param len length of seq.
 * @param k k-mer size.
 * @param large_wind_kmer_am
 * @param start 0-indexed position in seq to start hashing from. 
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
 * 
 * @throws BadThreadOutParams 
 */
// number of k-mers to be considered in the large window
template <class T>
void thread_wind(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
 */
template <class T>
void thread_wind(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * @param thread_count the number of threads to use
 * @param vec a vector of vectors in which the minimizers will be placed.
 *      Each vector corresponds to one thread. The minimizers within each vector
 *      will be in ascending order by index, and the vectors themselves will also
 *      be in ascending order by index, i.e. all minimizers in vector_i will go
 *      before all minimizers in vector_(i+1). 
 * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
 * @param len length of seq.
 * @param k k-mer size.
 * @param large_wind_kmer_am
 * @param start 0-indexed position in seq to start hashing from. 
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
 * 
 * @throws BadThreadOutParams 
 */
// number of k-mers to be considered in the large window
template <class T>
void thread_sync(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

/**
 * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
 */
template <class T>
void thread_sync(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

//------------- WORKER FUNCTIONS ----------------

// function that's passed to the thread for ModMinmizers
void thread_mod_roll(std::vector<uint32_t>& vec, const char* seq, 
    size_t ind, unsigned k, uint32_t mod, uint32_t congruence, 
    digest::MinimizedHashType minimized_h, unsigned assigned_kmer_am);

// function that's passed to the thread for WindowMinimizers
template <class T>
void thread_wind_roll(std::vector<uint32_t>& vec, const char* seq, 
    size_t ind, unsigned k, uint32_t large_wind_kmer_am, 
    digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am);

// function that's passed to the thread for Syncmers
template <class T>
void thread_sync_roll(std::vector<uint32_t>& vec, const char* seq, 
    size_t ind, unsigned k, uint32_t large_wind_kmer_am,
    digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am);

}


#include "thread_out.tpp"
#endif
