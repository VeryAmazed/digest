#include <array>
#include <benchmark/benchmark.h>
#include <cmath>
#include <cstdint>
#include <cwchar>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>
#include <set>
#include <deque>
#include <cassert>
#include <iostream>
#include <algorithm>

const int INPUT_SIZE = 1e5;

uint64_t inputs[INPUT_SIZE];
std::vector<int> st_outputs;
// std::vector<int> st_outputs2;
// std::vector<int> st_outputs3;
// std::vector<int> st_outputs4;
// std::vector<int> st_outputs5;
std::vector<int> st_outputs6;
// std::vector<int> mset_outputs;
// std::vector<int> naive_outputs;
// std::vector<int> naive_outputs2;
// std::vector<int> monoqueue_outputs;
// std::vector<int> monoqueue_outputs2;
// std::vector<int> monoqueue_outputs3;
// std::vector<int> monoqueue_outputs4;
// std::vector<int> monoqueue_outputs5;
// std::vector<int> monoqueue_outputs6;

void setupInput(){
    std::string path = "../tests/seg_tree/nums2.txt";
    std::ifstream ifs(path);
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    for(int i = 0; i < INPUT_SIZE; i++){
        ifs >> inputs[i];
    }
}

// class MonoQueue {
// 	uint64_t deq[64][2];
// 	int head = 0, tail = 0, i = 0;
// 	const int k;
//
// 	// cannot call this while full
// 	bool empty() {
// 		return head == tail;
// 	}
// public:
// 	MonoQueue(int k) : k(k) {}
//
// 	void add (uint64_t hash) {
// 		while (!empty() && deq[tail == 0 ? k - 1 : tail - 1][1] >= hash) {
// 			// pop back;
// 			if (--tail == -1) {
// 				tail = k - 1;
// 			}
// 		}
// 		if (!empty() && deq[head][0] == i - k) {
// 			// pop front
// 			if (++head == k) {
// 				head = 0;
// 			}
// 		}
//
// 		// append
// 		deq[tail][0] = i++;
// 		deq[tail][1] = hash;
// 		if (++tail == k) tail = 0;
// 	}
//
// 	int min() {
// 		return deq[head][0];
// 	}
// };







// class MonoQueue2 {
// 	__int128 deq[64];
//  	int head = 0, tail = 0, i = 0;
//  	const int k;
//
//  	// cannot call this while full
//  	bool empty() {
//  		return head == tail;
//  	}
//
// public:
//  	MonoQueue2(int k) : k(k) {}
//
//  	void add (uint64_t hash) {
//  		while (!empty() and deq[tail == 0 ? k - 1 : tail - 1] >> 32 >= hash) {
// 			if (--tail == -1) {
// 				tail = k - 1;
// 			}
//  		}
//  		if (!empty() and (int)deq[head] == i - k) {
// 			if (++head == k) {
// 				head = 0;
// 			}
//  		}
//
// 		deq[tail] = (__int128)hash << 32 | i++;
//  		if (++tail == k) tail = 0;
//  	}
//
//  	int min() {
//  		return deq[head];
//  	}
// };



// class MonoQueue3 {
// 	std::pair<uint64_t, int> deq[64];
//  	int head = 0, tail = 0, i = 0;
//  	const int k;
//
//  	// cannot call this while full
//  	bool empty() {
//  		return head == tail;
//  	}
//
// public:
//  	MonoQueue3(int k) : k(k) {}
//
//  	void add (uint64_t hash) {
//  		while (!empty() and deq[tail == 0 ? k - 1 : tail - 1].first >= hash) {
// 			if (--tail == -1) {
// 				tail = k - 1;
// 			}
//  		}
//  		if (!empty() and deq[head].second == i - k) {
// 			if (++head == k) {
// 				head = 0;
// 			}
//  		}
//
// 		deq[tail] = {hash, i++};
//  		if (++tail == k) tail = 0;
//  	}
//
//  	int min() {
//  		return deq[head].second;
//  	}
// };

// class MonoQueue4 {
// 	typedef struct {
// 		uint64_t hash;
// 		int index;
// 	} pair;
// 	pair deq[64];
//  	int head = 0, tail = 0, i = 0;
//  	const int k;
//
//  	// cannot call this while full
//  	bool empty() {
//  		return head == tail;
//  	}
//
// public:
//  	MonoQueue4(int k) : k(k) {}
//
//  	void add (uint64_t hash) {
//  		while (!empty() and deq[tail == 0 ? k - 1 : tail - 1].hash >= hash) {
// 			if (--tail == -1) {
// 				tail = k - 1;
// 			}
//  		}
//  		if (!empty() and deq[head].index == i - k) {
// 			if (++head == k) {
// 				head = 0;
// 			}
//  		}
//
// 		deq[tail] = {hash, i++};
//  		if (++tail == k) tail = 0;
//  	}
//
//  	int min() {
//  		return deq[head].index;
//  	}
// };

// class MonoQueue5 {
// 	std::vector<int> index;
// 	std::vector<uint64_t> hash;
//  	int head = 0, tail = 0, i = 0;
//  	const int k;
//
//  	// cannot call this while full
//  	bool empty() {
//  		return head == tail;
//  	}
//
// public:
//  	MonoQueue5(int k) : k(k) {
// 		index.reserve(k);
// 		hash.reserve(k);
// 	}
//
//  	void add(uint64_t hash) {
//  		while (!empty() and this->hash[tail == 0 ? k - 1 : tail - 1] >= hash) {
// 			if (--tail == -1) {
// 				tail = k - 1;
// 			}
//  		}
//  		if (!empty() and index[head] == i - k) {
// 			if (++head == k) {
// 				head = 0;
// 			}
//  		}
//
// 		this->hash[tail] = hash;
// 		index[tail] = i++;
//  		if (++tail == k) tail = 0;
//  	}
//
//  	int min() {
//  		return index[head];
//  	}
// };

// template<int k>
// struct MonoQueue6 {
//  	int head = 0, tail = 0;
// 	uint64_t queue[k];
//
// 	int timestamp[k];
// 	int time = 0;
//
//  	// cannot call this while full
//  	inline bool empty() {
//  		return head == tail;
//  	}
//
//  	void add (uint32_t index, uint32_t hash) {
// 		uint64_t q = (uint64_t)hash << 32 | index;
//  		while (!empty() and this->queue[tail == 0 ? k - 1 : tail - 1] >= q) {
// 			if (--tail == -1) {
// 				tail = k - 1;
// 			}
//  		}
//
//  		if (!empty() and this->timestamp[head] == time - k) {
// 			if (++head == k) {
// 				head = 0;
// 			}
//  		}
//
// 		this->queue[tail] = q;
// 		this->timestamp[tail] = time++;
//
//  		if (++tail == k) tail = 0;
//  	}
//
//  	uint32_t min() {
//  		return queue[head];
//  	}
// };


// static void BM_monoqueue_array(benchmark::State& state) {
//     for(auto _ : state) {
//         state.PauseTiming();
//         monoqueue_outputs.clear();
//         monoqueue_outputs.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
// 		MonoQueue mq(k);
//
// 		int i = 0;
// 		while (i < k) {
// 			mq.add(inputs[i++]);
// 		}
//
//         for (; i < INPUT_SIZE; i++){
// 			monoqueue_outputs.push_back(mq.min());
// 			mq.add(inputs[i]);
//         }
// 		monoqueue_outputs.push_back(mq.min());
// 	}
// }
// BENCHMARK(BM_monoqueue_array)->RangeMultiplier(2)->Range(1<<2, 1<<6);

// static void BM_monoqueue_128(benchmark::State& state) {
//     for(auto _ : state) {
//         state.PauseTiming();
//         monoqueue_outputs2.clear();
//         monoqueue_outputs2.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
// 		MonoQueue2 mq(k);
//
// 		int i = 0;
// 		while (i < k) {
// 			mq.add(inputs[i++]);
// 		}
//
//         for (; i < INPUT_SIZE; i++){
// 			monoqueue_outputs2.push_back(mq.min());
// 			mq.add(inputs[i]);
//         }
// 		monoqueue_outputs2.push_back(mq.min());
// 	}
// }
// BENCHMARK(BM_monoqueue_128)->RangeMultiplier(2)->Range(1<<2, 1<<6);

// static void BM_monoqueue_pair(benchmark::State& state) {
//     for(auto _ : state) {
//         state.PauseTiming();
//         monoqueue_outputs3.clear();
//         monoqueue_outputs3.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
// 		MonoQueue3 mq(k);
//
// 		int i = 0;
// 		while (i < k) {
// 			mq.add(inputs[i++]);
// 		}
//
//         for (; i < INPUT_SIZE; i++){
// 			monoqueue_outputs3.push_back(mq.min());
// 			mq.add(inputs[i]);
//         }
// 		monoqueue_outputs3.push_back(mq.min());
// 	}
// }
// BENCHMARK(BM_monoqueue_pair)->RangeMultiplier(2)->Range(1<<2, 1<<6);

// static void BM_monoqueue_struct(benchmark::State& state) {
//     for(auto _ : state) {
//         state.PauseTiming();
//         monoqueue_outputs4.clear();
//         monoqueue_outputs4.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
// 		MonoQueue4 mq(k);
//
// 		int i = 0;
// 		while (i < k) {
// 			mq.add(inputs[i++]);
// 		}
//
//         for (; i < INPUT_SIZE; i++){
// 			monoqueue_outputs4.push_back(mq.min());
// 			mq.add(inputs[i]);
//         }
// 		monoqueue_outputs4.push_back(mq.min());
// 	}
// }
// BENCHMARK(BM_monoqueue_struct)->RangeMultiplier(2)->Range(1<<2, 1<<6);


// static void BM_monoqueue_two_arrays(benchmark::State& state) {
//     for(auto _ : state) {
//         state.PauseTiming();
//         monoqueue_outputs5.clear();
//         monoqueue_outputs5.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
// 		MonoQueue5 mq(k);
//
// 		int i = 0;
// 		while (i < k) {
// 			mq.add(inputs[i++]);
// 		}
//
//         for (; i < INPUT_SIZE; i++){
// 			monoqueue_outputs5.push_back(mq.min());
// 			mq.add(inputs[i]);
//         }
// 		monoqueue_outputs5.push_back(mq.min());
// 	}
// }
// BENCHMARK(BM_monoqueue_two_arrays)->RangeMultiplier(2)->Range(1<<2, 1<<6);

// static void BM_monoqueue_template(benchmark::State& state) {
//     for(auto _ : state) {
//         state.PauseTiming();
//         monoqueue_outputs6.clear();
//         monoqueue_outputs6.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
/*
		 # define MONO() \
			int i = 0; \
			for (; i < k; i++) { \
				mq.add(i, inputs[i]); \
			} \
			for (; i < INPUT_SIZE; i++){ \
				monoqueue_outputs6.push_back(mq.min()); \
				mq.add(i, inputs[i]); \
			} \
			monoqueue_outputs6.push_back(mq.min());
*/
//
// 		switch (k) {
// 			case 2: {
// 				MonoQueue6<2> mq;
// 				MONO()
// 				break;
// 			}
// 			case 4: {
// 				MonoQueue6<4> mq;
// 				MONO()
// 				break;
// 			}
// 			case 8: {
// 				MonoQueue6<8> mq;
// 				MONO()
// 				break;
// 			}
// 			case 16: {
// 				MonoQueue6<16> mq;
// 				MONO()
// 				break;
// 			}
// 			case 32: {
// 				MonoQueue6<32> mq;
// 				MONO()
// 				break;
// 			}
// 			case 64: {
// 				MonoQueue6<64> mq;
// 				MONO()
// 				break;
// 			}
// 			default:
// 				throw std::runtime_error("monoqueue template");
// 		}
// 	}
// }
// BENCHMARK(BM_monoqueue_template)->RangeMultiplier(2)->Range(1<<2, 1<<6);












template <int k>
class MinSegmentTree {
	int i = 0;
	uint64_t segtree[2*k];

	constexpr int log2() {
		return std::ceil(std::log2(k));
	}
public:
	void set(uint32_t hash, int index) {
		int ind = i + k;
		if (++i == k) i = 0;

		segtree[ind] = (uint64_t)hash << 32 | index;
		for (int rep = log2(); rep >= 0; rep--) {
			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t min() {
		return segtree[1];
	}
};


static void BM_SegTree_template_unroll(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
        st_outputs.clear();
        st_outputs.reserve(INPUT_SIZE);
		const int k = state.range(0);

		// THIS IS A MACRO LOOK AWAY
		# define SEG() \
			int ind = 0; \
			while(ind < k){ \
				st.set(inputs[ind], ind); \
				ind++; \
			} \
			while(ind < INPUT_SIZE){ \
				st_outputs.push_back(st.min()); \
				st.set(inputs[ind], ind); \
				ind++; \
			} \
			st_outputs.push_back(st.min());

        state.ResumeTiming();

		if (k == 4) {
			MinSegmentTree<4> st;
			SEG()
		} else if (k == 8) {
			MinSegmentTree<8> st;
			SEG()
		} else if (k == 16) {
			MinSegmentTree<16> st;
			SEG()
		} else if (k == 32) {
			MinSegmentTree<32> st;
			SEG()
		} else if (k == 64) {
			MinSegmentTree<64> st;
			SEG()
		} else if (k == 128) {
			MinSegmentTree<128> st;
			SEG()
		} else {
			throw std::runtime_error("undeclared k-mer seg tree");
		}

		benchmark::DoNotOptimize(st_outputs);
		benchmark::DoNotOptimize(st_outputs.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_SegTree_template_unroll)->RangeMultiplier(2)->Range(1<<2, 1<<7);

template <int k>
class ArraySegTree {
	int i = 0;
	std::array<uint64_t,2*k> segtree;

	constexpr int log2() {
		return std::ceil(std::log2(k));
	}
public:
	void set(uint32_t hash, int index) {
		int ind = i + k;
		if (++i == k) i = 0;

		segtree[ind] = (uint64_t)hash << 32 | index;
		for (int rep = log2(); rep >= 0; rep--) {
			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t min() {
		return segtree[1];
	}
};


static void BM_SegTree_array(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
        st_outputs6.clear();
        st_outputs6.reserve(INPUT_SIZE);
		const int k = state.range(0);

		// THIS IS A MACRO LOOK AWAY
		# define SEG() \
			int ind = 0; \
			while(ind < k){ \
				st.set(inputs[ind], ind); \
				ind++; \
			} \
			while(ind < INPUT_SIZE){ \
				st_outputs6.push_back(st.min()); \
				st.set(inputs[ind], ind); \
				ind++; \
			} \
			st_outputs6.push_back(st.min());

        state.ResumeTiming();

		if (k == 4) {
			ArraySegTree<4> st;
			SEG()
		} else if (k == 8) {
			ArraySegTree<8> st;
			SEG()
		} else if (k == 16) {
			ArraySegTree<16> st;
			SEG()
		} else if (k == 32) {
			ArraySegTree<32> st;
			SEG()
		} else if (k == 64) {
			ArraySegTree<64> st;
			SEG()
		} else if (k == 128) {
			ArraySegTree<128> st;
			SEG()
		} else {
			throw std::runtime_error("undeclared k-mer seg tree");
		}

		benchmark::DoNotOptimize(st_outputs);
		benchmark::DoNotOptimize(st_outputs.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_SegTree_array)->RangeMultiplier(2)->Range(1<<2, 1<<7);

// template <int k>
// class MinSegmentTree0 {
// 	int i = 0;
// 	uint64_t segtree[2*k];
// public:
// 	void set(uint32_t hash, int index) {
// 		int ind = i + k;
// 		if (++i == k) i = 0;
//
// 		segtree[ind] = (uint64_t)hash << 32 | index;
// 		for (; ind > 1; ) {
// 			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
// 			ind >>= 1;
// 		}
// 	}
//
// 	uint32_t min() {
// 		return segtree[1];
// 	}
// };
//
// static void BM_SegTree_template(benchmark::State& state){
//     for(auto _ : state) {
//         state.PauseTiming();
//         st_outputs5.clear();
//         st_outputs5.reserve(INPUT_SIZE);
// 		const int k = state.range(0);
//
// 		// THIS IS A MACRO LOOK AWAY
/*
		# define SEG() \
			int ind = 0; \
			while(ind < k){ \
				st.set(inputs[ind], ind); \
				ind++; \
			} \
			while(ind < INPUT_SIZE){ \
				st_outputs5.push_back(st.min()); \
				st.set(inputs[ind], ind); \
				ind++; \
			} \
			st_outputs5.push_back(st.min());
*/
//
//         state.ResumeTiming();
//
// 		if (k == 4) {
// 			MinSegmentTree<4> st;
// 			SEG()
// 		} else if (k == 8) {
// 			MinSegmentTree<8> st;
// 			SEG()
// 		} else if (k == 16) {
// 			MinSegmentTree<16> st;
// 			SEG()
// 		} else if (k == 32) {
// 			MinSegmentTree<32> st;
// 			SEG()
// 		} else if (k == 64) {
// 			MinSegmentTree<64> st;
// 			SEG()
// 		} else if (k == 128) {
// 			MinSegmentTree<128> st;
// 			SEG()
// 		} else {
// 			throw std::runtime_error("undeclared k-mer seg tree");
// 		}
// 	}
// }
// BENCHMARK(BM_SegTree_template)->RangeMultiplier(2)->Range(1<<2, 1<<7);

// class MinSegmentTree2 {
// 	int i = 0;
// 	const int k, d;
// 	__int128 segtree[128];
//
// public:
// 	MinSegmentTree2(int k) : k(k), d(std::ceil(std::log2(k))) {}
//
// 	void set(uint64_t hash, int index) {
// 		int ind = i + k;
// 		if (++i == k) i = 0;
//
// 		segtree[ind] = (__int128)hash << 32 | index;
// 		for (int rep = d; rep >= 0; rep--) {
// 			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
// 			ind >>= 1;
// 		}
// 	}
//
// 	int min() {
// 		return segtree[1];
// 	}
// };
//
// static void BM_SegTree_128_unroll(benchmark::State& state){ // doesn't really unroll
//     for(auto _ : state) {
//         state.PauseTiming();
//         st_outputs2.clear();
//         st_outputs2.reserve(INPUT_SIZE);
// 		const int k = state.range(0);
//         state.ResumeTiming();
//
// 		MinSegmentTree2 st(k);
//
// 		int ind = 0;
// 		while(ind < k){
// 			st.set(inputs[ind], ind);
// 			ind++;
// 		}
// 		while(ind < INPUT_SIZE){
// 			st_outputs2.push_back(st.min());
// 			st.set(inputs[ind], ind);
// 			ind++;
// 		}
// 		st_outputs2.push_back(st.min());
// 	}
// }
// BENCHMARK(BM_SegTree_128_unroll)->RangeMultiplier(2)->Range(1<<2, 1<<6);


// class MinSegmentTree4 {
// 	int i = 0;
// 	const int k;
// 	uint64_t segtree[128];
// public:
// 	MinSegmentTree4(int k) : k(k) {}
//
// 	void set(uint32_t left, int right) {
// 		int ind = i + k;
// 		if (++i == k) i = 0;
//
// 		segtree[ind] = (uint64_t)left << 32 | right ;
// 		for (; ind > 1; ind >>= 1) {
// 			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
// 		}
// 	}
//
// 	uint32_t min() {
// 		return segtree[1];
// 	}
// };



// static void BM_SegTree_128(benchmark::State& state){
//     for(auto _ : state) {
//         state.PauseTiming();
//         st_outputs4.clear();
//         st_outputs4.reserve(INPUT_SIZE);
// 		const int k = state.range(0);
//         state.ResumeTiming();
//
//         int ind = 0;
//         MinSegmentTree4 st(k);
//         while(ind < k){
//             st.set(inputs[ind], ind);
//             ind++;
//         }
//         while(ind < INPUT_SIZE){
//             st_outputs4.push_back(st.min());
//             st.set(inputs[ind], ind);
//             ind++;
//         }
// 		st_outputs4.push_back(st.min());
// 	}
// }
// BENCHMARK(BM_SegTree_128)->RangeMultiplier(2)->Range(1<<2, 1<<6);

// class MinSegmentTree3 {
// 	int i = 0;
// 	const int k;
// 	std::pair<uint64_t,int> segtree[128];
// public:
// 	MinSegmentTree3(int k) : k(k) {}
//
// 	void set(uint64_t left, int right) {
// 		int ind = i + k;
// 		if (++i == k) i = 0;
//
// 		segtree[ind] = {left, right};
// 		for (; ind > 1; ind >>= 1) {
// 			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
// 		}
// 	}
//
// 	int min() {
// 		return segtree[1].second;
// 	}
// };
//
//
//
// static void BM_SegTree_pair(benchmark::State& state){
//     for(auto _ : state) {
//         state.PauseTiming();
//         st_outputs3.clear();
//         st_outputs3.reserve(INPUT_SIZE);
// 		const int k = state.range(0);
//         state.ResumeTiming();
//
//         int ind = 0;
//         MinSegmentTree3 st(k);
//         while(ind < k){
//             st.set(inputs[ind], ind);
//             ind++;
//         }
//         while(ind < INPUT_SIZE){
//             st_outputs3.push_back(st.min());
//             st.set(inputs[ind], ind);
//             ind++;
//         }
// 		st_outputs3.push_back(st.min());
// 	}
// }
// BENCHMARK(BM_SegTree_pair)->RangeMultiplier(2)->Range(1<<2, 1<<6);



















// static void BM_mset(benchmark::State& state){
//     for(auto _ : state) {
//         state.PauseTiming();
//         mset_outputs.clear();
//         mset_outputs.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         int ind = 0;
//         std::multiset<std::pair<int,int>> mset;
//         std::deque<decltype(mset)::iterator> vec;
//         while(ind < k){
//             vec.push_back(mset.emplace(inputs[ind], ind));
//             ind++;
//         }
//         int kick = 0;
//         state.ResumeTiming();
//         while(ind < INPUT_SIZE){
//             mset_outputs.push_back(mset.begin()->second);
//             mset.erase(vec[kick]);
//             vec[kick] = mset.emplace(inputs[ind], ind);
//             kick += 1;
//             if(kick == k) kick = 0;
//             ind++;
//         }
// 		mset_outputs.push_back(mset.begin()->second);
// 	}
// }
//BENCHMARK(BM_mset)->RangeMultiplier(2)->Range(1<<2, 1<<6);

// static void BM_naive(benchmark::State& state){
//     for(auto _ : state) {
//         state.PauseTiming();
//         naive_outputs.clear();
//         naive_outputs.reserve(INPUT_SIZE);
// 		int k = state.range(0);
//         state.ResumeTiming();
//
//         for(int i = 0; i+k <= INPUT_SIZE; i++){
//             uint64_t minAm = inputs[i];
// 			int ind = i;
//             for(int j = i+1; j < i + k; j++){
// 				if (inputs[j] < minAm) {
// 					minAm = inputs[j];
// 					ind = j;
// 				}
//             }
//             naive_outputs.push_back(ind);
//         }
// 	}
// }
// BENCHMARK(BM_naive)->RangeMultiplier(2)->Range(1<<2, 1<<6);


// static void BM_naive2(benchmark::State& state){
//     for(auto _ : state) {
//         state.PauseTiming();
//         naive_outputs2.clear();
//         naive_outputs2.reserve(INPUT_SIZE);
// 		const int k = state.range(0);
//         state.ResumeTiming();
//
// 		int cur = 0;
// 		for (int i = 0; i < k-1; i++) {
// 			if (inputs[i] < inputs[cur]) {
// 				cur = i;
// 			}
// 		}
//         for(int i = 0; i+k <= INPUT_SIZE; i++){
// 			if (inputs[i+k-1] < inputs[cur]) {
// 				cur = i+k-1;
// 			}
// 			if (cur == i-1) {
// 				cur++;
// 				for(int j = cur+1; j < i + k; j++){
// 					if (inputs[j] < inputs[cur]) {
// 						cur = j;
// 					}
// 				}
// 			}
//             naive_outputs2.push_back(cur);
//         }
// 	}
// }
// BENCHMARK(BM_naive2)->RangeMultiplier(2)->Range(1<<2, 1<<6);










int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    // sanity check
    // if it's correct for 1 window size, it's probably correct for all
    //assert(st_outputs == mset_outputs);
    // assert(st_outputs == naive_outputs);
    // assert(st_outputs == naive_outputs2);
    // assert(st_outputs == monoqueue_outputs);
    // assert(st_outputs == st_outputs2);
    // assert(st_outputs == st_outputs3);
    // assert(st_outputs == st_outputs4);
    // assert(st_outputs == st_outputs5);
    assert(st_outputs == st_outputs6);
    // assert(st_outputs == monoqueue_outputs2);
    // assert(st_outputs == monoqueue_outputs3);
    // assert(st_outputs == monoqueue_outputs4);
    // assert(st_outputs == monoqueue_outputs5);
    // assert(st_outputs == monoqueue_outputs6);
    std::cout << "Passed Asserts!" << std::endl;
}
