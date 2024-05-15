#ifndef WINDOW_MINIMIZER_HPP
#define WINDOW_MINIMIZER_HPP

#include "data_structure.hpp"
#include "digest/digester.hpp"
#include <cstddef>
#include <cstdint>

namespace digest {

/**
 * @brief Exception thrown when initializing a Window Minimizer or Syncmer with
 * a large window size of 0.
 *
 *
 */
class BadWindowSizeException : public std::exception {
	const char *what() const throw() {
		return "Number of kmers in large window cannot be 0";
	}
};

/**
 * @brief Child class of Digester that defines a minimizer as a kmer whose hash
 * is minimal among those in the large window. Parameters without a description
 * are the same as the parameters in the Digester parent class. They are simply
 * passed up to the parent constructor.
 *
 * @tparam P
 * @tparam T The data structure to use for performing range minimum queries to
 * find the minimal hash value.
 */
template <BadCharPolicy P, class T> class WindowMin : public Digester<P> {
  public:
	/**
	 * @param seq
	 * @param len
	 * @param k
	 * @param large_window the number of kmers in the large window, i.e. the
	 * number of kmers to be considered during the range minimum query.
	 * @param start
	 * @param minimized_h
	 *
	 * @throws BadWindowException thrown when large_window is passed in as 0
	 */
	WindowMin(const char *seq, size_t len, unsigned k, unsigned large_window,
			  size_t start = 0,
			  MinimizedHashType minimized_h = MinimizedHashType::CANON)
		: Digester<P>(seq, len, k, start, minimized_h), ds(large_window),
		  large_window(large_window), ds_size(0), is_minimized(false) {
		if (large_window == 0) {
			throw BadWindowSizeException();
		}
	}

	/**
	 * @param seq
	 * @param k
	 * @param large_window the number of kmers in the large window, i.e. the
	 * number of kmers to be considered during the range minimum query.
	 * @param start
	 * @param minimized_h
	 *
	 * @throws BadWindowException thrown when large_window is passed in as 0
	 */
	WindowMin(const std::string &seq, unsigned k, unsigned large_window,
			  size_t start = 0,
			  MinimizedHashType minimized_h = MinimizedHashType::CANON)
		: WindowMin<P, T>(seq.c_str(), seq.size(), k, large_window, start,
						  minimized_h) {}

	/**
	 * @brief adds up to amount of positions of minimizers into vec. Here a
	 * k-mer is considered a minimizer if its hash is the smallest in the large
	 * window. Rightmost index wins in ties
	 *
	 * @param amount
	 * @param vec
	 */
    void roll_minimizer(unsigned amount, std::vector<uint32_t>& vec) override {
		amount += vec.size();

		while (ds_size + 1 < large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				ds.insert(this->get_pos(), this->fhash);
			}else{
				ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			ds_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            roll_ds_wind(vec);
        }
    }

	/**
	 * @brief adds up to amount of positions and hashes of minimizers into vec.
	 * Here a k-mer is considered a minimizer if its hash is the smallest in the
	 * large window. Rightmost index wins in ties
	 *
	 * @param amount
	 * @param vec
	 */
	virtual void roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec) override {
		amount += vec.size();

		while (ds_size + 1 < large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				ds.insert(this->get_pos(), this->fhash);
			}else{
				ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			ds_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            roll_ds_wind(vec);
        }
    }



	/**
	 *
	 * @return unsigned, the value of large_window
	 */
	unsigned get_large_wind_kmer_am() { return large_window; }

	// function is mainly to help with tests
	size_t get_ds_size() { return ds_size; }

	// function is mainly to help with tests
	bool get_is_minimized() { return is_minimized; }

  protected:
	// data structure which will find miminum
	T ds;

	uint32_t large_window;

	// internal counter that tracks the number of actual values in the data
	// structure
	size_t ds_size;

	// internal bool keeping track of if we have obtained the first minimizer
	// yet, because we don't want to add a position to the vector if it's
	// already in there
	bool is_minimized;

	// the index of previous minimizer, a minimizer is only a new minimizer if
	// it is different from the previous minimizer
	uint32_t prev_mini;

  private:
	/**
	 * @brief helper function which handles adding the next hash into the data
	 * structure
	 *
	 */
    void roll_ds_wind(std::vector<uint32_t>& vec){
		if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
			ds.insert(this->get_pos(), this->chash);
		}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
			ds.insert(this->get_pos(), this->fhash);
		}else{
			ds.insert(this->get_pos(), this->rhash);
		}
		check(vec);
		
		this->roll_one();
    }

	/**
	 * @brief helper function which handles adding the next hash into the data
	 * structure
	 *
	 */
    void roll_ds_wind(std::vector<std::pair<uint32_t, uint32_t>>& vec){
		if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
			ds.insert(this->get_pos(), this->chash);
		}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
			ds.insert(this->get_pos(), this->fhash);
		}else{
			ds.insert(this->get_pos(), this->rhash);
		}
		check(vec);
		
		this->roll_one();
    }


	/**
	 * @brief helper function that checks to see if the current minimizer is a
	 * new minimizer, and should thus be added to the vec
	 *
	 * @param vec
	 */
    void check(std::vector<uint32_t>& vec){
        if(is_minimized){
            if(ds.min() != prev_mini){
                prev_mini =	ds.min();
                vec.emplace_back(prev_mini);
            }
        }else{
            is_minimized = true;
            prev_mini = ds.min();
			vec.emplace_back(prev_mini);
        }
    }


	/**
	 * @brief helper function that checks to see if the current minimizer is a
	 * new minimizer, and should thus be added to the vec
	 *
	 * @param vec
	 */
    void check(std::vector<std::pair<uint32_t, uint32_t>>& vec){
        if(is_minimized){
            if(ds.min() != prev_mini){
                prev_mini =	ds.min();
                vec.emplace_back(prev_mini, ds.min_hash());
            }
        }else{
            is_minimized = true;
            prev_mini = ds.min();
			vec.emplace_back(prev_mini, ds.min_hash());
        }
    }
};

} // namespace digest

#endif // WINDOW_MINIMIZER_HPP
