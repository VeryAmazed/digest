#ifndef MOD_MINI_HPP
#define MOD_MINI_HPP

#include "digest/digester.hpp"
#include <cassert>
#include <cstdint>

namespace digest {

/**
 * @brief Exception thrown when initializing a mod minimizer object where the
 * target value after modding is greater than the mod value.
 *
 *
 */
class BadModException : public std::exception {
	const char *what() const throw() {
		return "mod must be greater than congruence.";
	}
};

/**
 * @brief Child class of Digester that defines a minimizer as a kmer whose hash
 * is equal to some target value after being modded. Parameters without a
 * description are the same as the parameters in the Digester parent class. They
 * are simply passed up to the parent constructor.
 *
 * @tparam P
 */
template <BadCharPolicy P> class ModMin : public Digester<P> {
  public:
	/**
	 * @brief
	 *
	 * @param seq
	 * @param len
	 * @param k
	 * @param mod mod space to be used to calculate universal minimizers
	 * @param congruence value we want minimizer hashes to be congruent to in
	 * the mod space
	 * @param start
	 * @param minimized_h
	 *
	 * @throws BadModException Thrown when congruence is greater or equal to mod
	 */
	ModMin(const char *seq, size_t len, unsigned k, uint32_t mod,
		   uint32_t congruence = 0, size_t start = 0,
		   MinimizedHashType minimized_h = MinimizedHashType::CANON)
		: Digester<P>(seq, len, k, start, minimized_h), mod(mod),
		  congruence(congruence) {
		if (congruence >= mod) {
			throw BadModException();
		}
	}

	/**
	 *
	 * @param seq
	 * @param k
	 * @param mod mod space to be used to calculate universal minimizers
	 * @param congruence value we want minimizer hashes to be congruent to in
	 * the mod space
	 * @param start
	 * @param minimized_h
	 *
	 * @throws BadModException Thrown when congruence is greater or equal to mod
	 */
	ModMin(const std::string &seq, unsigned k, uint32_t mod,
		   uint32_t congruence = 0, size_t start = 0,
		   MinimizedHashType minimized_h = MinimizedHashType::CANON)
		: ModMin<P>(seq.c_str(), seq.size(), k, mod, congruence, start,
					minimized_h) {}

	/**
	 * @brief adds up to amount of positions of minimizers into vec. Here a
	 * k-mer is considered a minimizer if its hash is congruent to congruence in
	 * the mod space.
	 *
	 * @param amount
	 * @param vec
	 */
	void roll_minimizer(unsigned amount, std::vector<uint32_t> &vec) override;

	/**
	 * @brief adds up to amount of positions and hashes of minimizers into vec.
	 * Here a k-mer is considered a minimizer if its hash is congruent to
	 * congruence in the mod space.
	 *
	 * @param amount
	 * @param vec
	 */
	void
	roll_minimizer(unsigned amount,
				   std::vector<std::pair<uint32_t, uint32_t>> &vec) override;

	/**
	 * @return uint32_t, the mod space being used
	 */
	uint32_t get_mod() { return mod; }

	/**
	 * @return uint32_t, the value the minimized hash must be congruent to
	 */
	uint32_t get_congruence() { return congruence; }

  private:
	uint32_t mod;
	uint32_t congruence;
};

} // namespace digest

#include "mod_minimizer.tpp"
#endif
