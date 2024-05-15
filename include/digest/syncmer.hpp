#ifndef SYNCMER_HPP
#define SYNCMER_HPP

#include "digest/digester.hpp"
#include "digest/window_minimizer.hpp"

namespace digest {

/**
 * @brief This class inherits from WindowMinimizer (implementation reasons), but
 * the represent very different things. A Syncmer is defined as a large window
 * where the minimal hash among all kmers in the large window belong to either
 * the leftmost or rightmost kmer. Parameters without a description are the same
 * as the parameters in the Digester parent class. They are simply passed up to
 * the parent constructor.
 *
 * @tparam P
 * @tparam T The data structure to use for performing range minimum queries to
 * find the minimal hash value.
 */
template <BadCharPolicy P, class T> class Syncmer : public WindowMin<P, T> {
  public:
	/**
	 *
	 * @param seq
	 * @param len
	 * @param k
	 * @param large_window the number of kmers in the large window, i.e. the
	 * number of kmers to be considered during the range minimum query.
	 * @param start
	 * @param minimized_h
	 *
	 * @throws BadWindowException Thrown when large_window is passed in as 0
	 */
	Syncmer(const char *seq, size_t len, unsigned k, unsigned large_window,
			size_t start = 0,
			MinimizedHashType minimized_h = MinimizedHashType::CANON)
		: WindowMin<P, T>(seq, len, k, large_window, start, minimized_h) {}

	/**
	 *
	 * @param seq
	 * @param k
	 * @param large_window the number of kmers in the large window, i.e. the
	 * number of kmers to be considered during the range minimum query.
	 * @param start
	 * @param minimized_h
	 *
	 * @throws BadWindowException Thrown when large_window is passed in as 0
	 */
	Syncmer(const std::string &seq, unsigned k, unsigned large_window,
			size_t start = 0,
			MinimizedHashType minimized_h = MinimizedHashType::CANON)
		: Syncmer<P, T>(seq.c_str(), seq.size(), k, large_window, start,
						minimized_h) {}

	/**
	 * @brief adds up to amount of positions of syncmers into vec. Here
	 * a large window is considered a syncmer if the smallest hash in the large
	 * window is at the leftmost or rightmost position.
	 *
	 * @param amount
	 * @param vec
	 */
	void roll_minimizer(unsigned amount, std::vector<uint32_t> &vec) {
		amount += vec.size();

		while (this->ds_size + 1 < this->large_window and this->is_valid_hash) {
			if (this->get_minimized_h() == digest::MinimizedHashType::CANON) {
				this->ds.insert(this->get_pos(), this->chash);
			} else if (this->get_minimized_h() ==
					   digest::MinimizedHashType::FORWARD) {
				this->ds.insert(this->get_pos(), this->fhash);
			} else {
				this->ds.insert(this->get_pos(), this->rhash);
			}

			this->roll_one();
			this->ds_size++;
		}

		while (this->is_valid_hash and vec.size() < amount) {
			Syncmer::roll_ds_sync(vec);
		}
	}

	/**
	 * @brief adds up to amount of positions and hashes of syncmers into vec.
	 * Here a large window is considered a syncmer if the smallest hash in the
	 * large window is at the leftmost or rightmost position.
	 *
	 * @param amount
	 * @param vec
	 */
	void roll_minimizer(unsigned amount,
						std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		amount += vec.size();

		while (this->ds_size + 1 < this->large_window and this->is_valid_hash) {
			if (this->get_minimized_h() == digest::MinimizedHashType::CANON) {
				this->ds.insert(this->get_pos(), this->chash);
			} else if (this->get_minimized_h() ==
					   digest::MinimizedHashType::FORWARD) {
				this->ds.insert(this->get_pos(), this->fhash);
			} else {
				this->ds.insert(this->get_pos(), this->rhash);
			}

			this->roll_one();
			this->ds_size++;
		}

		while (this->is_valid_hash and vec.size() < amount) {
			Syncmer::roll_ds_sync(vec);
		}
	}

  private:
	/**
	 * @brief helper function which handles adding the next hash into the data
	 * structure
	 *
	 */
	void roll_ds_sync(std::vector<uint32_t> &vec) {
		if (this->get_minimized_h() == digest::MinimizedHashType::CANON) {
			this->ds.insert(this->get_pos(), this->chash);
		} else if (this->get_minimized_h() ==
				   digest::MinimizedHashType::FORWARD) {
			this->ds.insert(this->get_pos(), this->fhash);
		} else {
			this->ds.insert(this->get_pos(), this->rhash);
		}
		this->ds.min_syncmer(vec);

		this->roll_one();
	}

	/**
	 * @brief helper function which handles adding the next hash into the data
	 * structure
	 *
	 */
	void roll_ds_sync(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		if (this->get_minimized_h() == digest::MinimizedHashType::CANON) {
			this->ds.insert(this->get_pos(), this->chash);
		} else if (this->get_minimized_h() ==
				   digest::MinimizedHashType::FORWARD) {
			this->ds.insert(this->get_pos(), this->fhash);
		} else {
			this->ds.insert(this->get_pos(), this->rhash);
		}
		this->ds.min_syncmer(vec);

		this->roll_one();
	}
};

} // namespace digest

#endif // SYNCMER_HPP
