#ifndef THREAD_OUT_HPP
#define THREAD_OUT_HPP

#include "digest/mod_minimizer.hpp"
#include "digest/window_minimizer.hpp"
#include "digest/syncmer.hpp"
#include <thread>

/*
    Possible implementation for multi-threading the digestion of a single sequence.
    The key thing to note is basically by carefully telling where each digester should start digesting
    you can make it so each kmer is only considered once.
    I have very little experience with threading, so you could probably thread this out better than me
*/
namespace thread_out{

class BadThreadOutParams : public std::exception
{
	const char * what () const throw ()
    {
    	return "k must be greater than 3, start must be less than len, \
        and num threads must be greater or equal to the number of kmers";
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
 * @throws BadThreadOutParams Thrown if k is less than 4, start is after the end of the string, 
 *      or the number of threads exceeds the number of possible kmers in the string
 */
void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t start = 0, 
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
 * @param mod mod space to be used to calculate universal minimizers
 * @param congruence value we want minimizer hashes to be congruent to in the mod space
 * @param start 0-indexed position in seq to start hashing from. 
 * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
 * 
 * @throws BadThreadOutParams Thrown if k is less than 4, start is after the end of the string, 
 *      or the number of threads exceeds the number of possible kmers in the string
 */
void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const std::string& seq, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

void thread_wind(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const char* seq, size_t len, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

void thread_wind(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const std::string& seq, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

void thread_sync(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const char* seq, size_t len, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

void thread_sync(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const std::string& seq, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

// function that's passed to the thread for ModMinmizers
void thread_mod_roll(std::vector<size_t>& vec, const char* seq, 
    size_t ind, unsigned k, uint64_t mod, uint64_t congruence, 
    digest::MinimizedHashType minimized_h, unsigned assigned_kmer_am);

}

#endif