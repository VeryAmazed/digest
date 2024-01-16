#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <sys/types.h>
#include <vector>
#include <utility>
#include <stdint.h>

// requirement on all data_structures
// constructor which accepts uint32_t
// void set(uint32_t index, (uint 32/64) hash)
// uint32_t min() // returns minimum
// pair<uint32_t, uint 32/64> min_hash() // returns index, hash
// void min_syncmer(vector<uint 32/64> &vec) // appends minimum if syncmer
// void min_syncmer(vector<pair<uint 32, uint 32>> &vec) // appends (left syncmer index, right syncmer index) if syncmer
// assignment/copy constructors if you want to use them

namespace data_structure {

// template<uint32_t k>
// struct MonoQueue {
//  	int head = 0, tail = 0;
// 	// one extra slot so empty() works when full
// 	// {hash, index, time}
// 	std::array<std::array<uint32_t,k+1>,3> queue;
// 	uint32_t time = 0;
//
//  	bool empty() {
//  		return head == tail;
//  	}
//
// 	MonoQueue(uint32_t) {
// 		queue[2].fill(0); // necessary to avoid 1 in a billion collision
// 		queue[0].fill(0); // gets rid of warning
// 	}
// 	MonoQueue(const MonoQueue& other) = default;
// 	MonoQueue &operator=(const MonoQueue& other) = default;
//
//  	void insert(uint32_t index, uint32_t hash) {
//  		if (queue[2][head] == time - k) {
// 			if (++head == k+1) {
// 				head = 0;
// 			}
//  		}
//
//  		while (not empty() and queue[0][tail == 0 ? k : tail - 1] >= hash) {
// 			if (--tail == -1) {
// 				tail = k;
// 			}
//  		}
//
// 		queue[0][tail] = hash;
// 		queue[1][tail] = index;
// 		queue[2][tail] = time++;
//
//  		if (++tail == k+1) tail = 0;
//
//  	}
//
// 	uint32_t min() {
//  		return queue[1][head];
// 	}
//
// 	void min_syncmer(std::vector<uint32_t> &vec) {
// 		if (queue[2][head] == time - k or queue[2][head] == time - 1) {
// 			vec.emplace_back(queue[1][head]);
// 		}
// 	}
// };


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

	void insert(uint32_t index, uint32_t hash) {
		int ind = i + k;
		if (++i == k) i = 0;

		// negate so we can use max so that ties are broken by rightmost
		segtree[ind] = (uint64_t)~hash << 32 | index;
		for (int rep = 0; rep < log2(); rep++) {
			segtree[ind >> 1] = std::max(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t min() {
		return segtree[1];
	}

	std::pair<uint32_t, uint32_t> min_hash() {
		return {segtree[1], ~(segtree[1] >> 32)};
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (segtree[1] == segtree[i] or segtree[1] == segtree[i ? i-1 : k-1]) {
			vec.emplace_back(segtree[i]);
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t,uint32_t>> &vec) {
		if (segtree[1] == segtree[i] or segtree[1] == segtree[i ? i-1 : k-1]) {
			vec.emplace_back(segtree[i], segtree[i ? i-1 : k-1]);
		}
	}
};

// template<uint32_t k>
// struct Set {
// 	std::set<uint64_t> mset;
// 	std::array<std::set<uint64_t>::iterator,k> vec;
// 	int i = 0;
//
// 	Set(uint32_t) {
// 		// edge case where hash == all 1's can cause a set collision.
// 		for (unsigned i = 0; i < k; i++) {
// 			vec[i] = mset.emplace(i).first;
// 		}
// 	}
// 	Set(const Set& other) = delete; // have to copy over iterators
// 	Set &operator=(const Set& other) = delete;
//
// 	void insert(uint32_t index, uint32_t hash) {
// 		mset.erase(vec[i]);
//
// 		vec[i] = mset.emplace((uint64_t)~hash << 32 | index).first;
// 		if (++i == k) i = 0;
// 	}
//
// 	uint32_t min() {
// 		return *mset.rbegin();
// 	}
//
// 	void min_syncmer(std::vector<uint32_t> &vec) {
//
// 	}
// };

template<uint32_t k>
struct Naive {
	std::array<uint64_t,k> arr;
	uint i = 0;

	Naive(uint32_t) {};
	Naive(const Naive& other) = default;
	Naive &operator=(const Naive& other) = default;

	void insert(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)~hash << 32 | index;
		if (++i == k) i = 0;
	}

	uint32_t min() {
		int i = k-1;
		for (int j = k-2; j >= 0; j--) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return arr[i];
	}

	std::pair<uint32_t, uint32_t> min_hash() {
		int i = k-1;
		for (int j = k-2; j >= 0; j--) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return {arr[i], ~(uint32_t)(arr[i] >> 32)};
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		uint j = k-1;
		for (int l = k-2; l >= 0; l--) {
			if (arr[l] > arr[j]) {
				j = l;
			}
		}
		if (j == i or j == (i ? i-1 : k-1)) {
			vec.emplace_back(arr[i]);
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t,uint32_t>> &vec) {
		uint j = k-1;
		for (int l = k-2; l >= 0; l--) {
			if (arr[l] > arr[j]) {
				j = l;
			}
		}
		if (j == i or j == (i ? i-1 : k-1)) {
			vec.emplace_back(arr[i], arr[i ? i-1 : k-1]);
		}
	}
};

template<uint32_t k>
struct Naive2 {
	uint i = 0;
	uint last = 0;
	std::vector<uint64_t> arr = std::vector<uint64_t>(k);

	Naive2(uint32_t) {};
	Naive2(const Naive2& other) = default;
	Naive2 &operator=(const Naive2& other) = default;

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

		if (++i == k) i = 0;
	}

	uint32_t min() {
		return arr[last];
	}

	std::pair<uint32_t, uint32_t> min_hash() {
		return {arr[last], ~(uint32_t)(arr[last] >> 32)};
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (last == i or last == (i ? i-1 : k-1)) {
			vec.emplace_back(arr[i]);
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t,uint32_t>> &vec) {
		if (last == i or last == (i ? i-1 : k-1)) {
			vec.emplace_back(arr[i], arr[i ? i-1 : k-1]);
		}
	}
};

struct Adaptive {
	uint32_t k, i = 0, last = 0;
	std::vector<uint64_t> arr;

	Adaptive(uint32_t k) : k(k), arr(k) {}
	Adaptive(const Adaptive& other) = default;
	Adaptive &operator=(const Adaptive& other) = default;

	void naive(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)~hash << 32 | index;
		if (++i == k) i = 0;
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

		if (++i == k) i = 0;
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
			int i = k-1;
			for (int j = k-2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return arr[i];
		} else {
			return arr[last];
		}
	}

	std::pair<uint32_t, uint32_t> min_hash() {
		if (k < 16) {
			int i = k-1;
			for (int j = k-2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return {arr[i], ~(uint32_t)(arr[i] >> 32)};
		} else {
			return {arr[last], ~(uint32_t)(arr[last] >> 32)};
		}
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (k < 16) {
			uint j = k-1;
			for (int l = k-2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (j == i or j == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i]);
			}
		} else {
			if (last == i or last == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i]);
			}
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		if (k < 16) {
			uint j = k-1;
			for (int l = k-2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (j == i or j == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i], arr[i ? i-1 : k-1]);
			}
		} else {
			if (last == i or last == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i], arr[i ? i-1 : k-1]);
			}
		}
	}
};

struct Adaptive64 {
	uint32_t k, i = 0, last = 0;
	std::vector<__uint128_t> arr;

	Adaptive64(uint32_t k) : k(k), arr(k) {}
	Adaptive64(const Adaptive64& other) = default;
	Adaptive64 &operator=(const Adaptive64& other) = default;

	void naive(uint32_t index, uint64_t hash) {
		arr[i] = (__uint128_t)~hash << 32 | index;
		if (++i == k) i = 0;
	}

	void naive2(uint32_t index, uint64_t hash) {
		// flip the hash bits so we can take the maximum
		arr[i] = (__uint128_t)~hash << 32 | index;

		if (arr[i] > arr[last]) {
			last = i;
		} else if (last == i) {
			for (int j = k-1; j >= 0; j--) {
				if (arr[j] > arr[last]) {
					last = j;
				}
			}
		}

		if (++i == k) i = 0;
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
			int i = k-1;
			for (int j = k-2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return arr[i];
		} else {
			return arr[last];
		}
	}

	std::pair<uint32_t, uint64_t> min_hash() {
		if (k < 16) {
			int i = k-1;
			for (int j = k-2; j >= 0; j--) {
				if (arr[j] > arr[i]) {
					i = j;
				}
			}
			return {arr[i], ~(uint64_t)(arr[i] >> 32)};
		} else {
			return {arr[last], ~(uint64_t)(arr[last] >> 32)};
		}
	}

	void min_syncmer(std::vector<uint32_t> &vec) {
		if (k < 16) {
			uint j = k-1;
			for (int l = k-2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (j == i or j == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i]);
			}
		} else {
			if (last == i or last == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i]);
			}
		}
	}

	void min_syncmer(std::vector<std::pair<uint32_t, uint32_t>> &vec) {
		if (k < 16) {
			uint j = k-1;
			for (int l = k-2; l >= 0; l--) {
				if (arr[l] > arr[j]) {
					j = l;
				}
			}
			if (j == i or j == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i], arr[i ? i-1 : k-1]);
			}
		} else {
			if (last == i or last == (i ? i-1 : k-1)) {
				vec.emplace_back(arr[i], arr[i ? i-1 : k-1]);
			}
		}
	}
};

}
