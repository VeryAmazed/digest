#pragma once

#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>
#include <stdint.h>


namespace segtree{
// Based on a template taken from USACO.guide and then modified by me (for competitive programming), and now modified again (for this)
// https://usaco.guide/gold/PURS?lang=cpp
// https://codeforces.com/blog/entry/18051 (USACO.guide was probably heavily inspired by this)
/** A data structure that can answer point update & range minimum queries. */

// k is size of segtree
template <uint32_t k>
struct SegTree {
	int i = 0;
	std::array<uint64_t,2*k> segtree;

	// ceil(log2(k))
	constexpr int log2() {
		assert(k > 0);
		int ans = 31 - __builtin_clz(k-1); // floor(log2(k-1))
		return ans + 1;
	}

	void set(uint32_t hash, uint32_t index) {
		int ind = i + k; 
		if (++i == k) i = 0;

		// we store the xor of the index so larger indices are favored
		segtree[ind] = (uint64_t)hash << 32 | (0xffffffff ^ index);
		for (int rep = 0; rep < log2(); rep++) {
			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	// ties are won by rightmost index
	uint32_t min() {
		return segtree[1] ^ 0xffffffff;
	}

	uint32_t get_hash(uint32_t actual_index) {
		return segtree[actual_index] >> 32;
	}
	int32_t get_index(uint32_t actual_index) {
		return segtree[actual_index] ^ 0xffffffff;
	}
};

}
