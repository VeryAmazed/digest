#include <array>
#include <benchmark/benchmark.h>
#include <cstdint>
#include <random>
#include <iostream>
#include <limits.h>
#include <sys/types.h>
#include <vector>

const int INPUT_SIZE = 1e7;
std::array<uint32_t,2*INPUT_SIZE> hashes;

std::map<int,std::array<uint32_t,INPUT_SIZE>> st_outputs;
std::map<int,std::array<uint32_t,INPUT_SIZE>> st_outputs_normal;
std::map<int,std::array<uint32_t,INPUT_SIZE>> mset_outputs;
std::map<int,std::array<uint32_t,INPUT_SIZE>> naive_outputs;
std::map<int,std::array<uint32_t,INPUT_SIZE>> naive_outputs2;
std::map<int,std::array<uint32_t,INPUT_SIZE>> monoqueue_outputs;

void setupInput(){
	std::random_device rd; // seed
	std::mt19937 gen(rd()); // generator
	std::uniform_int_distribution<uint32_t> distrib(0,UINT_MAX); // [0, 2**32]
	for (uint32_t &h : hashes) {
        h = distrib(gen);
    }
}

template<uint32_t k>
class MonoQueue {
 	int head = 0, tail = 0;
	// one extra slot so empty() works when full
	// {hash, index, time}
	std::array<std::array<uint32_t,k+1>,3> queue;
	uint32_t time = 0;

 	bool empty() {
 		return head == tail;
 	}

public:
	MonoQueue() {
		queue[2].fill(0); // necessary to avoid 1 in a billion collision
		queue[0].fill(0); // gets rid of warning
	}

 	void add (uint32_t index, uint32_t hash) {
 		if (queue[2][head] == time - k) {
			if (++head == k+1) {
				head = 0;
			}
 		}

 		while (not empty() and queue[0][tail == 0 ? k : tail - 1] <= hash) {
			if (--tail == -1) {
				tail = k;
			}
 		}

		queue[0][tail] = hash;
		queue[1][tail] = index;
		queue[2][tail] = time++;

 		if (++tail == k+1) tail = 0;
 	}

 	uint32_t max() {
 		return queue[1][head];
 	}
};

template<int k>
static void BM_monoqueue(benchmark::State& state) {
	auto &out = monoqueue_outputs[k];
    for(auto _ : state) {
		MonoQueue<k> mq;
		for (int i = 0; i < k; i++) {
			mq.add(i, hashes[i]);
		}
		for (int i = 0; i < INPUT_SIZE; i++){
			out[i] = mq.max();
			mq.add(i+k, hashes[i+k]);
		}
	}
}

template<int K=1024>
class SegmentTree {
	int i = 0;
	const int k;
	uint64_t segtree[2*K];

	constexpr int log2() {
		return std::ceil(std::log2(K));
	}
public:
	SegmentTree(int k) : k(k) {}
	void set(uint32_t hash, uint32_t index) {
		int ind = i + k;
		if (++i == k) i = 0;

		segtree[ind] = (uint64_t)hash << 32 | index;
		for (int rep = 0; rep < log2(); rep++) {
			segtree[ind >> 1] = std::max(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t max() {
		return segtree[1];
	}
};

template<int k>
static void BM_segtree(benchmark::State& state) {
	auto &out = st_outputs[k];
    for(auto _ : state) {
		SegmentTree<k> st(k);

		for (int i = 0; i < k; i++){
			st.set(hashes[i], i);
		}
		for (int i = 0; i < INPUT_SIZE; i++){
			out[i] = st.max();
			st.set(hashes[i+k], i+k);
		}
	}
}

class Segtree_normal {
	int i = 0;
	const int k;
	std::vector<uint64_t> segtree;

	// constexpr int log2() {
	// 	return std::ceil(std::log2(k));
	// }
public:
	Segtree_normal(int k) : k(k), segtree(2*k) {}

	void set(uint32_t hash, uint32_t index) {
		int ind = i + k;
		if (++i == k) i = 0;

		segtree[ind] = (uint64_t)hash << 32 | index;
		while (ind > 1) {
			segtree[ind >> 1] = std::max(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t max() {
		return segtree[1];
	}
};

template<int k>
static void BM_segtree_normal(benchmark::State& state) {
	auto &out = st_outputs_normal[k];
    for(auto _ : state) {
		Segtree_normal st(k);

		for (int i = 0; i < k; i++){
			st.set(hashes[i], i);
		}
		for (int i = 0; i < INPUT_SIZE; i++){
			out[i] = st.max();
			st.set(hashes[i+k], i+k);
		}
	}
}

template<uint32_t k>
static void BM_mset(benchmark::State& state){
	auto &out = mset_outputs[k];
    for(auto _ : state) {
        std::multiset<std::uint64_t> mset;
        std::array<decltype(mset)::iterator,k> vec;
		for (uint32_t i = 0; i < k; i++) {
            vec[i] = mset.emplace((uint64_t)hashes[i] << 32 | i);
        }

        for (uint32_t i = 0, kick = 0; i < INPUT_SIZE; i++){
            out[i] = *mset.rbegin();

            mset.erase(vec[kick]);
            vec[kick] = mset.emplace((uint64_t)hashes[i+k] << 32 | (i+k));
            if (++kick == k) kick = 0;
        }
	}
}

class Naive {
	std::vector<uint64_t> arr;
	int i = 0;
	const int k;

public:
	Naive(int k) : arr(k), k(k) {}

	void insert(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)hash << 32 | index;
		if (++i == k) i = 0;
	}

	uint32_t max() {
		int i = 0;
		for (int j = 1; j < k; j++) {
			if (arr[j] > arr[i]) {
				i = j;
			}
		}
		return arr[i];
	}
};

template<int k>
static void BM_naive(benchmark::State& state){
	auto &out = naive_outputs[k];
    for(auto _ : state) {
		Naive n(k);
		for (int i = 0; i < k; i++) {
			n.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            out[i] = n.max();
			n.insert(i+k, hashes[i+k]);
        }
	}
}

class Naive2 {
	const int k;
	int i = 0;
	int last = 0;
	std::vector<uint64_t> arr;

public:
	Naive2(int k) : k(k), arr(k) {}

	void insert(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)hash << 32 | index;

		if (arr[i] > arr[last]) {
			last = i;
		} else if (last == i) {
			for (int j = 0; j < k; j++) {
				if (arr[j] > arr[last]) {
					last = j;
				}
			}
		}

		if (++i == k) i = 0;
	}

	uint32_t max() {
		return arr[last];
	}
};

template<int k>
static void BM_naive2(benchmark::State& state){
	auto &out = naive_outputs2[k];
    for(auto _ : state) {
		Naive2 n(k);
		for (int i = 0; i < k; i++) {
			n.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            out[i] = n.max();
			n.insert(i+k, hashes[i+k]);
        }
	}
}

// template<int k>
// class Naive2 {
// 	std::array<std::array<uint32_t, k>,2> arr; // hash, index
// 	int i = 0;
// 	int last = 0;
//
// public:
//
// 	void insert(uint32_t index, uint32_t hash) {
// 		arr[0][i] = hash;
// 		arr[1][i] = index;
//
// 		if (arr[0][i] > arr[0][last]) {
// 			last = i;
// 		} else if (last == i) {
// 			for (int j = i-1; j >= 0; j--) {
// 				if (arr[0][j] > arr[0][last]) {
// 					last = j;
// 				}
// 			}
// 			for (int j = k-1; j > i; j--) {
// 				if (arr[0][j] > arr[0][last]) {
// 					last = j;
// 				}
// 			}
// 		}
//
// 		if (++i == k) i = 0;
// 	}
//
// 	uint32_t max() {
// 		return arr[1][last];
// 	}
// };
//
// template<int k>
// static void BM_naive2(benchmark::State& state){
// 	auto &out = naive_outputs2[k];
//     for(auto _ : state) {
// 		Naive2<k> n;
// 		for (int i = 0; i < k; i++) {
// 			n.insert(i, hashes[i]);
// 		}
//         for (int i = 0; i < INPUT_SIZE; i++) {
//             out[i] = n.max();
// 			n.insert(i+k, hashes[i+k]);
//         }
// 	}
// }

#define test(name) \
	BENCHMARK_TEMPLATE(name, 2); \
	BENCHMARK_TEMPLATE(name, 3); \
	BENCHMARK_TEMPLATE(name, 4); \
	BENCHMARK_TEMPLATE(name, 6); \
	BENCHMARK_TEMPLATE(name, 8); \
	BENCHMARK_TEMPLATE(name, 12); \
	BENCHMARK_TEMPLATE(name, 16); \
	BENCHMARK_TEMPLATE(name, 24); \
	BENCHMARK_TEMPLATE(name, 32); \
	BENCHMARK_TEMPLATE(name, 48); \
	BENCHMARK_TEMPLATE(name, 64); \
	BENCHMARK_TEMPLATE(name, 96); \
	BENCHMARK_TEMPLATE(name, 128); \
	BENCHMARK_TEMPLATE(name, 256); \
	BENCHMARK_TEMPLATE(name, 512); \
	BENCHMARK_TEMPLATE(name, 1024); \

test(BM_segtree);
test(BM_segtree_normal);
test(BM_naive);
test(BM_naive2);

// test(BM_monoqueue);
// test(BM_mset);

int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    // sanity check
    assert(st_outputs == naive_outputs);
    assert(st_outputs == naive_outputs2);
    assert(st_outputs == st_outputs_normal);
    // assert(st_outputs == monoqueue_outputs);
    // assert(st_outputs == mset_outputs);
    std::cout << "Passed Asserts!" << std::endl;
}
