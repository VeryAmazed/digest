
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

class WindowMin : public Digester{
	public:
		/**
         * @param seq 
         * @param len
         * @param k 
         * @param large_wind_kmer_am number of kmers to be considered at a time
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        WindowMin(const char* seq, size_t len, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, unsigned minimized_h = 0)
        :  Digester(seq, len, k, start, minimized_h), st(segtree::SegTree(large_wind_kmer_am)), large_wind_kmer_am(large_wind_kmer_am), st_index(0), st_size(0), is_minimized(false)
        {	
            if(large_wind_kmer_am == 0){
				throw BadWindowSizeException();
			}
        }

        /**
         * @param seq 
         * @param k 
         * @param large_wind_kmer_am number of kmers to be considered at a time
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        WindowMin(const std::string& seq, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, unsigned minimized_h = 0) :
            WindowMin(seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h)
        {}

		WindowMin(const WindowMin& copy) : Digester(copy), st(segtree::SegTree(copy.st)), large_wind_kmer_am(copy.large_wind_kmer_am), st_index(copy.st_index), st_size(copy.st_size), is_minimized(copy.is_minimized), prev_mini(copy.prev_mini)
		{
		}

		WindowMin& operator=(const WindowMin& copy){
			this->large_wind_kmer_am = copy.large_wind_kmer_am;
			this->st_index = copy.st_index;
			this->st_size = copy.st_size;
			this->is_minimized = copy.is_minimized;
			this->prev_mini = copy.prev_mini;
			st = copy.st;
			Digester::operator=(copy);
			return *this;
		}

		virtual ~WindowMin(){
		}

		/**
		 * @brief adds up to amount of positions of minimizers into vec, here a k-mer is considered a minimizer if its hash is the smallest in the large window, using rightmost index wins in ties 
         *        Time Complexity: O(log(large_wind_kmer_am)) per k-mer tested
		 * 
		 * @param amount 
		 * @param vec 
		 */
		virtual void roll_minimizer(unsigned amount, std::vector<size_t>& vec) override;

		unsigned get_large_wind_kmer_am(){
			return large_wind_kmer_am;
		}
		
		size_t get_st_index(){
			return st_index;
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
		segtree::SegTree st;

		// number of k-mer to be considered in the large window
		unsigned large_wind_kmer_am;

		// internal index that denotes the position of the leftmost element in the large window within the segment tree
		size_t st_index;

		// internal counter that tracks the number of actual values in the segment tree
		size_t st_size;

		// internal bool keeping track of if we have obtained the first minimizer yet, because we don't want to add a position to the vector if it's already in there
		bool is_minimized;

		// internal pair representing the previous minimizer, a minimizer is only a new minimizer if it is different from the previous minimizer
		std::pair<uint64_t, size_t> prev_mini;

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
#endif
