
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
	// extra stuff needed, segtree, large window size, segtree internal pointer (private)
	// errors, large window size can't be 0
	public:
		/**
         * 
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
        :  Digester(seq, len, k, start, minimized_h), large_wind_kmer_am(large_wind_kmer_am), st_index(0), st_size(0), is_minimized(false)
        {	
            if(large_wind_kmer_am == 0){
				throw BadWindowSizeException();
			}
			(this->st) = new segtree::SegTree(large_wind_kmer_am);
        }

        /**
         * 
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

		WindowMin(const WindowMin& copy) : Digester(copy), large_wind_kmer_am(copy.large_wind_kmer_am), st_index(copy.st_index), st_size(copy.st_size), is_minimized(copy.is_minimized), prev_mini(copy.prev_mini)
		{
			st = new segtree::SegTree(*copy.st);
		}

		WindowMin& operator=(const WindowMin& copy){
			this->large_wind_kmer_am = copy.large_wind_kmer_am;
			this->st_index = copy.st_index;
			this->st_size = copy.st_size;
			this->is_minimized = copy.is_minimized;
			this->prev_mini = copy.prev_mini;
			*st = *(copy.st);
			Digester::operator=(copy);
			return *this;
		}

		virtual ~WindowMin(){
			delete st;
		}

		virtual void roll_minimizer(unsigned amount, std::vector<size_t>& vec) override;

		unsigned get_large_wind_kmer_am(){
			return large_wind_kmer_am;
		}

		// I don't know why you would ever call this, but it's here
		size_t get_st_index(){
			return st_index;
		}

		size_t get_st_size(){
			return st_size;
		}

		bool get_is_minimized(){
			return is_minimized;
		}

	protected:
		segtree::SegTree* st;
		unsigned large_wind_kmer_am;
		size_t st_index;
		size_t st_size;
		bool is_minimized;
		std::pair<uint64_t, size_t> prev_mini;

		void fill_st();
	
	private:
		void check(std::vector<size_t>& vec);
};

}
#endif
