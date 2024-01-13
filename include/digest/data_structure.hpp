#pragma once

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <utility>
#include <stdint.h>

// requirement on all data_structures
// constructor which accepts uint32_t
// uint32_t set(uint32_t index, uint32_t hash) // returns new minimum
// assignment/copy constructors if you want to use them

namespace data_structure {

template<uint32_t k>
struct MonoQueue {
 	int head = 0, tail = 0;
	// one extra slot so empty() works when full
	// {hash, index, time}
	std::array<std::array<uint32_t,k+1>,3> queue;
	uint32_t time = 0;

 	bool empty() {
 		return head == tail;
 	}

	MonoQueue(uint32_t) {
		queue[2].fill(0); // necessary to avoid 1 in a billion collision
		queue[0].fill(0); // gets rid of warning
	}
	MonoQueue(const MonoQueue& other) = default;
	MonoQueue &operator=(const MonoQueue& other) = default;

 	uint32_t insert(uint32_t index, uint32_t hash) {
 		if (queue[2][head] == time - k) {
			if (++head == k+1) {
				head = 0;
			}
 		}

 		while (not empty() and queue[0][tail == 0 ? k : tail - 1] >= hash) {
			if (--tail == -1) {
				tail = k;
			}
 		}

		queue[0][tail] = hash;
		queue[1][tail] = index;
		queue[2][tail] = time++;

 		if (++tail == k+1) tail = 0;

 		return queue[1][head];
 	}
};


// Based on a template taken from USACO.guide and then modified by me (for competitive programming), and now modified again (for this)
// https://usaco.guide/gold/PURS?lang=cpp
// https://codeforces.com/blog/entry/18051 (USACO.guide was probably heavily inspired by this)
/** A data structure that can answer point update & range minimum queries. */
template<int k>
struct SegmentTree {
	int i = 0;
	std::array<uint64_t,2*k> segtree;

	constexpr int log2() {
		return std::ceil(std::log2(k));
	}

	SegmentTree(uint32_t) {}
	SegmentTree(const SegmentTree& other) = default;
	SegmentTree &operator=(const SegmentTree& other) = default;

	uint32_t insert(uint32_t index, uint32_t hash) {
		int ind = i + k;
		if (++i == k) i = 0;

		// negate so we can use max so that ties are broken by rightmost
		segtree[ind] = (uint64_t)~hash << 32 | index;
		for (int rep = 0; rep < log2(); rep++) {
			segtree[ind >> 1] = std::max(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}

		return segtree[1];
	}
};

template<uint32_t k>
struct Set {
	std::set<uint64_t> mset;
	std::array<std::set<uint64_t>::iterator,k> vec;
	int i = 0;

	Set(uint32_t) {
		// edge case where hash == all 1's can cause a set collision.
		for (unsigned i = 0; i < k; i++) {
			vec[i] = mset.emplace(i).first;
		}
	}
	Set(const Set& other) = delete; // have to copy over iterators
	Set &operator=(const Set& other) = delete;

	uint32_t insert(uint32_t index, uint32_t hash) {
		mset.erase(vec[i]);

		vec[i] = mset.emplace((uint64_t)~hash << 32 | index).first;
		if (++i == k) i = 0;

		return *mset.rbegin();
	}
};

template<uint32_t k>
struct Naive {
	std::array<uint64_t,k> arr;
	int i = 0;

	Naive(uint32_t) {};
	Naive(const Naive& other) = default;
	Naive &operator=(const Naive& other) = default;

	uint32_t insert(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)~hash << 32 | index;
		if (++i == k) i = 0;

		// get min
		int i = k-1;
		for (int j = k-2; j >= 0; j--) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return arr[i];
	}
};

template<uint32_t k>
struct Naive2 {
	int i = 0;
	int last = 0;
	std::vector<uint64_t> arr = std::vector<uint64_t>(k);

	Naive2(uint32_t) {};
	Naive2(const Naive2& other) = default;
	Naive2 &operator=(const Naive2& other) = default;

	uint32_t insert(uint32_t index, uint32_t hash) {
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

		if (++i == k) i = 0;

		return arr[last];
	}
};

struct Adaptive {
	uint32_t k, i = 0, last = 0;
	std::vector<uint64_t> arr;

	Adaptive(uint32_t k) : k(k), arr(k) {}
	Adaptive(const Adaptive& other) = default;
	Adaptive &operator=(const Adaptive& other) = default;

	uint32_t naive(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)~hash << 32 | index;
		if (++i == k) i = 0;

		// get min
		int i = k-1;
		for (int j = k-2; j >= 0; j--) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return arr[i];
	}

	uint32_t naive2(uint32_t index, uint32_t hash) {
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

		if (++i == k) i = 0;

		return arr[last];
	}

	uint32_t insert(uint32_t index, uint32_t hash) {
		if (k < 16) {
			return naive(index, hash);
		}
		return naive2(index, hash);
	}
};

}
