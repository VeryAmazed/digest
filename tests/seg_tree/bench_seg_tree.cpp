#include <array>
#include <benchmark/benchmark.h>
#include <cstdint>
#include <random>
#include <iostream>
#include <limits.h>
#include <vector>

const int INPUT_SIZE = 1e7;
std::array<uint32_t,2*INPUT_SIZE> hashes;

std::map<int,std::array<uint32_t,INPUT_SIZE>> st_outputs;
std::map<int,std::array<uint32_t,INPUT_SIZE>> mset_outputs;
std::map<int,std::array<uint32_t,INPUT_SIZE>> naive_outputs;
std::map<int,std::array<uint32_t,INPUT_SIZE>> naive_outputs1;
std::map<int,std::array<uint32_t,INPUT_SIZE>> naive_outputs2;
std::map<int,std::array<uint32_t,INPUT_SIZE>> naive_outputs2c;
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

template<int k>
class SegmentTree {
	int i = 0;
	uint64_t segtree[2*k];

	constexpr int log2() {
		return std::ceil(std::log2(k));
	}
public:
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
static void BM_Segment_Tree(benchmark::State& state) {
	auto &out = st_outputs[k];
    for(auto _ : state) {
		SegmentTree<k> st;

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

template<int k>
class Naive {
	std::array<uint64_t, k> arr;
	int i = 0;

public:
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
		Naive<k> n;
		for (int i = 0; i < k; i++) {
			n.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            out[i] = n.max();
			n.insert(i+k, hashes[i+k]);
        }
	}
}

template<int k>
class Naive1 {
	std::array<std::array<uint32_t, k>,2> arr; // hash, index
	int i = 0;

public:
	void insert(uint32_t index, uint32_t hash) {
		arr[0][i] = hash;
		arr[1][i] = index;
		if (++i == k) i = 0;
	}

	uint32_t max() {
		int i = k-1;
		for (int j = k-2; j >= 0; j--) {
			if (arr[0][j] > arr[0][i]) {
				i = j;
			}
		}
		return arr[1][i];
	}
};

template<int k>
static void BM_naive1(benchmark::State& state){
	auto &out = naive_outputs1[k];
    for(auto _ : state) {
		Naive1<k> n;
		for (int i = 0; i < k; i++) {
			n.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            out[i] = n.max();
			n.insert(i+k, hashes[i+k]);
        }
	}
}

template<int k>
class Naive2Combined {
	std::array<uint64_t, k> arr;
	int i = 0;
	int last = 0;

public:
	Naive2Combined() {
		arr.fill(0); // gets rid of warning
	}

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
static void BM_naive2_combined(benchmark::State& state){
	auto &out = naive_outputs2c[k];
    for(auto _ : state) {
		Naive2Combined<k> n;
		for (int i = 0; i < k; i++) {
			n.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            out[i] = n.max();
			n.insert(i+k, hashes[i+k]);
        }
	}
}

template<int k>
class Naive2 {
	std::array<uint64_t, k> arr; // hash, index
	int i = 0;
	int last = 0;

public:

	void insert(uint32_t index, uint32_t hash) {
		arr[i] = (uint64_t)hash << 32 | index; // this is faster

		if (arr[i]>>32 > arr[last]>>32) {
			last = i;
		} else if (last == i) {
			for (int j = i-1; j >= 0; j--) {
				if (arr[j]>>32 > arr[last]>>32) {
					last = j;
				}
			}
			for (int j = k-1; j > i; j--) {
				if (arr[j]>>32 > arr[last]>>32) {
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
		Naive2<k> n;
		for (int i = 0; i < k; i++) {
			n.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            out[i] = n.max();
			n.insert(i+k, hashes[i+k]);
        }
	}
}

#define test(name) \
	BENCHMARK_TEMPLATE(name, 1); \
	BENCHMARK_TEMPLATE(name, 2); \
	BENCHMARK_TEMPLATE(name, 4); \
	BENCHMARK_TEMPLATE(name, 8); \
	BENCHMARK_TEMPLATE(name, 16); \
	BENCHMARK_TEMPLATE(name, 32); \
	BENCHMARK_TEMPLATE(name, 64); \
	BENCHMARK_TEMPLATE(name, 128); \
	BENCHMARK_TEMPLATE(name, 256); \
	BENCHMARK_TEMPLATE(name, 512); \
	BENCHMARK_TEMPLATE(name, 1024); \

test(BM_monoqueue);
test(BM_Segment_Tree);
test(BM_naive);
test(BM_naive1);
test(BM_naive2);
test(BM_naive2_combined);
test(BM_mset);

int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    assert(naive_outputs2 == naive_outputs2c);
    assert(naive_outputs1 == naive_outputs);
    // sanity check
    assert(st_outputs == mset_outputs);
    assert(st_outputs == naive_outputs);
    assert(st_outputs == naive_outputs2);
    assert(st_outputs == monoqueue_outputs);
    std::cout << "Passed Asserts!" << std::endl;
}
