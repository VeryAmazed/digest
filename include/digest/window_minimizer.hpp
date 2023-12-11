
#ifndef WIND_MIN_HPP
#define WIND_MIN_HPP
#include "digester.hpp"
#include "point_update_st.hpp"

namespace digest{

class BadWindowSizeException : public std::exception
{
	const char * what () const throw ()
    {
    	return "Number of kmers in large window cannot be 0";
    }
};

// number of k-mer to be considered in the large window
template <int32_t large_window>
class WindowMin : public Digester{
	public:
		/**
         * @param seq 
         * @param len
         * @param k 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        WindowMin(const char* seq, size_t len, unsigned k, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON)
        :  Digester(seq, len, k, start, minimized_h), st_size(0), is_minimized(false)
        {	
            if(large_window == 0){
				throw BadWindowSizeException();
			}
        }

        /**
         * @param seq 
         * @param k 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        WindowMin(const std::string& seq, unsigned k, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON) :
            WindowMin(seq.c_str(), seq.size(), k, start, minimized_h)
        {}

		/**
		 * @brief adds up to amount of positions of minimizers into vec, here a k-mer is considered a minimizer if its hash is the smallest in the large window, using rightmost index wins in ties 
         *        Time Complexity: O(log(large_wind_kmer_am)) per k-mer tested
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
		// internal point update segment tree data structure used to find minimum in large window
		segtree::SegTree<large_window> st;

		// internal counter that tracks the number of actual values in the segment tree
		size_t st_size;

		// internal bool keeping track of if we have obtained the first minimizer yet, because we don't want to add a position to the vector if it's already in there
		bool is_minimized;

		// the index of previous minimizer, a minimizer is only a new minimizer if it is different from the previous minimizer
		int32_t prev_mini;

		/**
		 * @brief helper function which handles adding new elements into the segment tree when it is not full
		 * 
		 */
		void fill_st();
	
	private:

		/**
		 * @brief helper function that checks to see if the current minimizer is a new minimizer, and should thus be added to the vec
		 * 
		 * @param vec 
		 */
		void check(std::vector<size_t>& vec);
};

}

#include "window_minimizer.tpp"
#endif
