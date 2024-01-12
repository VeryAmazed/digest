#pragma once
#include "digester.hpp"
#include "data_structure.hpp"
#include <cstddef>
#include <cstdint>

namespace digest{

class BadWindowSizeException : public std::exception
{
	const char * what () const throw ()
    {
    	return "Number of kmers in large window cannot be 0";
    }
};

// number of k-mers to be considered in the large window
template <class T>
class WindowMin : public Digester{
	public:
		/**
         * @param seq 
         * @param len
         * @param k 
         * @param large_window 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        WindowMin(const char* seq, size_t len, unsigned k, unsigned large_window, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON)
        :  Digester(seq, len, k, start, minimized_h), ds(large_window), st_size(0), is_minimized(false)
        {	
            if(large_window == 0){
				throw BadWindowSizeException();
			}
        }

        /**
         * @param seq 
         * @param k 
         * @param large_window 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        WindowMin(const std::string& seq, unsigned k, unsigned large_window, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON) :
            WindowMin(seq.c_str(), seq.size(), k, large_window, start, minimized_h)
        {}

		/**
		 * @brief adds up to amount of positions of minimizers into vec, here a k-mer is considered a minimizer if its hash is the smallest in the large window, using rightmost index wins in ties 
		 * 
		 * @param amount 
		 * @param vec 
		 */
		virtual void roll_minimizer(unsigned amount, std::vector<size_t>& vec) override;

		unsigned get_large_wind_kmer_am(){
			return large_window;
		}
		
		size_t get_st_size(){
			return st_size;
		}

		// function is mainly to help with tests
		bool get_is_minimized(){
			return is_minimized;
		}

	protected:
		// data structure which will find miminum
		T ds;

		uint32_t large_window;

		// internal counter that tracks the number of actual values in the segment tree
		size_t st_size;

		// internal bool keeping track of if we have obtained the first minimizer yet, because we don't want to add a position to the vector if it's already in there
		bool is_minimized;

		// the index of previous minimizer, a minimizer is only a new minimizer if it is different from the previous minimizer
		uint32_t prev_mini;

		/**
		 * @brief helper function which handles adding new elements into the segment tree when it is not full
		 * 
		 */
		void fill_st(std::vector<size_t>& vec);
	
	private:

		/**
		 * @brief helper function that checks to see if the current minimizer is a new minimizer, and should thus be added to the vec
		 * 
		 * @param vec 
		 */
		void check(std::vector<size_t>& vec, uint32_t hash);
};

}

#include "window_minimizer.tpp"
