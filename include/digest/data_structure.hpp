#ifndef DATA_STRUCTURE_HPP
#define DATA_STRUCTURE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdint.h>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * Data structures for minimum hash queries on a window.
 * ntHash does not support `large_window < 4`
 *
 * Selecting the correct `data_structure`
 * our general guidelines:
 * * for `large_window` < 12, use Naive
 * * for 12 <= `large_window` <= 16 use SegmentTree
 * * for `large_window` > 16 use Naive2
 *
 * Adaptive performs at worst about 10% slower than best
 * Adaptive64 performs at worst about 100% slower than best
 */

namespace digest::ds {

/**
 * All data_structures must follow this interface. Add the min_syncmer functions
 * for syncmer support.
 */
template <typename T> struct Interface {
	static_assert(std::is_same<T, uint32_t>() || std::is_same<T, uint64_t>(),
				  "T must be either uint32_t or uint64_t");

	/** constructor must accept uint32_t large_window */
	Interface(uint32_t);

	/** returns the index of the minimum hash */
	virtual uint32_t min();
	/** returns the minimum hash */
	virtual T min_hash();

	/** appends minimum if syncmer */
	virtual void min_syncmer(std::vector<uint32_t> &vec);
	/** appends (left syncmer index, right syncmer index) */
	virtual void min_syncmer(std::vector<std::pair<uint32_t, T>> &vec);
};

// Based on a template taken from USACO.guide and then modified by me (for
// competitive programming), and now modified again (for this)
// https://usaco.guide/gold/PURS?lang=cpp
// https://codeforces.com/blog/entry/18051 (USACO.guide was probably heavily
// inspired by this)
/** A data structure that can answer point update & range minimum queries. */

/**
 * @brief Segment Tree data structure. Supports log(n) point updates and range
 * minimum queries.
 *
 * @tparam k large window size
 */
template <int k> struct SegmentTree {
	int i = k;
	std::array<uint64_t, 2 * k> segtree = {};

	constexpr int log2() { return std::ceil(std::log2(k)); }

	SegmentTree(uint32_t) {}
	SegmentTree(const SegmentTree &other) = default;
	SegmentTree &operator=(const SegmentTree &other) = default;

	void insert(uint32_t index, uint32_t hash) {
		int ind = i;
		if (++i == 2 * k)
			i = k;

		// negate so we can use max so that ties are broken by rightmost
		segtree[ind] = (uint64_t)~hash << 32 | index;
		for (int rep = 0; rep < log2(); rep++) {
			segtree[ind >> 1] = std::max(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t min() { return segtree[1]; }

	uint32_t min_hash() { return ~(segtree[1] >> 32); }

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (segtree[1] >> 32 ==
			std::max(uint32_t(segtree[i] >> 32),
					 uint32_t(segtree[i == k ? 2 * k - 1 : i - 1] >> 32))) {
			vec.emplace_back(segtree[i]);
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		if (segtree[1] >> 32 ==
			std::max(uint32_t(segtree[i] >> 32),
					 uint32_t(segtree[i == k ? 2 * k - 1 : i - 1] >> 32))) {
			vec.emplace_back(segtree[i], ~(segtree[1] >> 32));
		}
	}
};

/**
 * @brief Naive data structure. Naively loops through the array to find the
 * minimum.
 *
 * @tparam k large window size
 */
template <uint32_t k> struct Naive {
	std::array<uint64_t, k> arr;
	unsigned int i = 0;

	Naive(uint32_t){};
	Naive(const Naive &other) = default;
	Naive &operator=(const Naive &other) = default;

	void insert(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)~hash << 32 | index;
		if (++i == k)
			i = 0;
	}

	uint32_t min() {
		int i = k - 1;
		for (int j = k - 2; j >= 0; j--) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return arr[i];
	}

	uint32_t min_hash() {
		int i = k - 1;
		for (int j = k - 2; j >= 0; j--) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return ~(uint32_t)(arr[i] >> 32);
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		unsigned int j = 0;
		for (unsigned int l = 1; l < k; l++) {
			if (arr[l] > arr[j]) {
				j = l;
			}
		}
		if (arr[j] >> 32 == std::max(uint32_t(arr[i] >> 32),
									 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
			vec.emplace_back(arr[i]);
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		unsigned int j = k - 1;
		for (int l = k - 2; l >= 0; l--) {
			if (arr[l] > arr[j]) {
				j = l;
			}
		}
		if (arr[j] >> 32 == std::max(uint32_t(arr[i] >> 32),
									 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
			vec.emplace_back(arr[i], ~(uint32_t)(arr[j] >> 32));
		}
	}
};

/**
 * @brief Naive2 data structure. Remembers the last minimum index and only loops
 * through the array when this index leaves the window.
 *
 * @tparam k large window size
 */
template <uint32_t k> struct Naive2 {
	unsigned int i = 0;
	unsigned int last = 0;
	std::vector<uint64_t> arr = std::vector<uint64_t>(k);

	Naive2(uint32_t) {};
	Naive2(const Naive2 &other) = default;
	Naive2 &operator=(const Naive2 &other) = default;

	void insert(uint32_t index, uint32_t hash) {
		// flip the hash bits so we can take the maximum
		arr[i] = (uint64_t)~hash << 32 | index;

		if (arr[i] > arr[last]) {
			last = i;
		} else if (last == i) {
			for (unsigned j = 0; j < k; j++) {
				if (arr[j] > arr[last]) {
					last = j;
				}
			}
		}

		if (++i == k)
			i = 0;
	}

	uint32_t min() { return arr[last]; }

	uint32_t min_hash() { return ~(uint32_t)(arr[last] >> 32); }

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (arr[last] >> 32 ==
			std::max(uint32_t(arr[i] >> 32),
					 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
			vec.emplace_back(arr[i]);
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		if (arr[last] >> 32 ==
			std::max(uint32_t(arr[i] >> 32),
					 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
			vec.emplace_back(arr[i], ~(uint32_t)(arr[last] >> 32));
		}
	}
};

/**
 * @brief Adaptive data structure. Selects between Naive and Naive2 based on the
 * large window size.
 */
struct Adaptive {
	uint32_t k, i = 0, last = 0;
	std::vector<uint64_t> arr;

	Adaptive(uint32_t k) : k(k), arr(k) {}
	Adaptive(const Adaptive &other) = default;
	Adaptive &operator=(const Adaptive &other) = default;

	void naive(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)~hash << 32 | index;
		if (++i == k)
			i = 0;
	}

	void naive2(uint32_t index, uint32_t hash) {
		// flip the hash bits so we can take the maximum
		arr[i] = (uint64_t)~hash << 32 | index;

		if (arr[i] > arr[last]) {
			last = i;
		} else if (last == i) {
			for (unsigned j = 0; j < k; j++) {
				if (arr[j] > arr[last]) {
					last = j;
				}
			}
		}

		if (++i == k)
			i = 0;
	}

	void insert(uint32_t index, uint32_t hash) {
		if (k < 16) {
			naive(index, hash);
		} else {
			naive2(index, hash);
		}
	}

	uint32_t min() {
		if (k < 16) {
			int i = k - 1;
			for (int j = k - 2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return arr[i];
		} else {
			return arr[last];
		}
	}

	uint32_t min_hash() {
		if (k < 16) {
			int i = k - 1;
			for (int j = k - 2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return ~(uint32_t)(arr[i] >> 32);
		} else {
			return ~(uint32_t)(arr[last] >> 32);
		}
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (k < 16) {
			unsigned int j = k - 1;
			for (int l = k - 2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (arr[j] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i]);
			}
		} else {
			if (arr[last] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i]);
			}
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		if (k < 16) {
			unsigned int j = k - 1;
			for (int l = k - 2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (arr[j] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i], ~(uint32_t)(arr[j] >> 32));
			}
		} else {
			if (arr[last] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i], ~(uint32_t)(arr[last] >> 32));
			}
		}
	}
};

/**
 * @brief Same as Adaptive but uses 64-bit hashes.
 */
struct Adaptive64 {
	uint32_t k, i = 0, last = 0;
	std::vector<__uint128_t> arr;

	Adaptive64(uint32_t k) : k(k), arr(k) {}
	Adaptive64(const Adaptive64 &other) = default;
	Adaptive64 &operator=(const Adaptive64 &other) = default;

	void naive(uint32_t index, uint64_t hash) {
		arr[i] = (__uint128_t)~hash << 32 | index;
		if (++i == k)
			i = 0;
	}

	void naive2(uint32_t index, uint64_t hash) {
		// flip the hash bits so we can take the maximum
		arr[i] = (__uint128_t)~hash << 32 | index;

		if (arr[i] > arr[last]) {
			last = i;
		} else if (last == i) {
			for (int j = k - 1; j >= 0; j--) {
				if (arr[j] > arr[last]) {
					last = j;
				}
			}
		}

		if (++i == k)
			i = 0;
	}

	void insert(uint32_t index, uint64_t hash) {
		if (k < 16) {
			naive(index, hash);
		} else {
			return naive2(index, hash);
		}
	}

	uint32_t min() {
		if (k < 16) {
			int i = k - 1;
			for (int j = k - 2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return arr[i];
		} else {
			return arr[last];
		}
	}

	uint64_t min_hash() {
		if (k < 16) {
			int i = k - 1;
			for (int j = k - 2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return ~(uint64_t)(arr[i] >> 32);
		} else {
			return ~(uint64_t)(arr[last] >> 32);
		}
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (k < 16) {
			unsigned int j = k - 1;
			for (int l = k - 2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (arr[j] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i]);
			}
		} else {
			if (arr[last] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i]);
			}
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint64_t>> &vec) {
		if (k < 16) {
			unsigned int j = k - 1;
			for (int l = k - 2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (arr[j] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i], ~(uint64_t)(arr[j] >> 32));
			}
		} else {
			if (arr[last] >> 32 ==
				std::max(uint32_t(arr[i] >> 32),
						 uint32_t(arr[i ? i - 1 : k - 1] >> 32))) {
				vec.emplace_back(arr[i], ~(uint64_t)(arr[last] >> 32));
			}
		}
	}
};
} // namespace digest::ds

#endif // DATA_STRUCTURE_HPP
