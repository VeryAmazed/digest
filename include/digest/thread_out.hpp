#ifndef THREAD_OUT_HPP
#define THREAD_OUT_HPP

#include "digest/mod_minimizer.hpp"
#include "digest/window_minimizer.hpp"
#include "digest/syncmer.hpp"
#include <thread>

namespace thread_out{

class BadThreadOutParams : public std::exception
{
	const char * what () const throw ()
    {
    	return "k must be greater than 3, start must be less than len, \
        and num threads must be greater or equal to the number of kmers";
    }
};

void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t start = 0, 
    digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON);

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

void worker_roll(digest::Digester& dig, std::vector<size_t>& vec, unsigned amount);

}

#endif