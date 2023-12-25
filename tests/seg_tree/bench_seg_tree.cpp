#include <bits/stdc++.h>
#include <benchmark/benchmark.h>

const int INPUT_SIZE = 1e6;
std::array<uint32_t,INPUT_SIZE> hashes;

std::array<uint32_t,INPUT_SIZE> st_outputs;
std::array<uint32_t,INPUT_SIZE> st_outputs2;
std::array<uint32_t,INPUT_SIZE> mset_outputs;
std::array<uint32_t,INPUT_SIZE> naive_outputs;
std::array<uint32_t,INPUT_SIZE> monoqueue_outputs6;

// nums seems to be upto 1e9 (length = 1e5)
// nums2 is 64-bit hashes (length = 1e5)
// nums3 is 32-bit hashes (length = 1e6)

void setupInput(){
    std::string path = "../tests/seg_tree/nums3.txt";
    std::ifstream ifs(path);
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    for(int i = 0; i < INPUT_SIZE; i++){
        ifs >> hashes[i];
    }
}

// TODO:
// test segtree with uint32_t[2] array, pair, struct

template<int k>
struct MonoQueue6 {
 	int head = 0, tail = 0;
	// one extra slot so empty() works when full
	std::array<uint64_t,k+1> queue;

	int timestamp[k+1];
	int time = 0;

 	inline bool empty() {
 		return head == tail;
 	}

 	void add (uint32_t index, uint32_t hash) {
		uint64_t q = (uint64_t)hash << 32 | ~index;
 		while (!empty() and this->queue[tail == 0 ? k : tail - 1] >= q) {
			if (--tail == -1) {
				tail = k;
			}
 		}

 		if (!empty() and this->timestamp[head] == time - k) {
			if (++head == k+1) {
				head = 0;
			}
 		}

		this->queue[tail] = q;
		this->timestamp[tail] = time++;

 		if (++tail == k+1) tail = 0;
 	}

 	uint32_t min() {
 		return ~(uint32_t)queue[head];
 	}
};

template<uint32_t k>
void MONO(MonoQueue6<k> &mq) {
	int i = 0;
	for (; i < k; i++) {
		mq.add(i, hashes[i]);
	}
	for (; i < INPUT_SIZE; i++){
		monoqueue_outputs6[i-k] = mq.min();
		mq.add(i, hashes[i]);
	}
	monoqueue_outputs6[i-k] = mq.min();
}

static void BM_monoqueue_template(benchmark::State& state) {
    for(auto _ : state) {
		switch (state.range(0)) {
			case 2: {
				MonoQueue6<2> mq;
				MONO<2>(mq);
				break;
			}
			case 4: {
				MonoQueue6<4> mq;
				MONO<4>(mq);
				break;
			}
			case 8: {
				MonoQueue6<8> mq;
				MONO<8>(mq);
				break;
			}
			case 16: {
				MonoQueue6<16> mq;
				MONO<16>(mq);
				break;
			}
			case 32: {
				MonoQueue6<32> mq;
				MONO<32>(mq);
				break;
			}
			case 64: {
				MonoQueue6<64> mq;
				MONO<64>(mq);
				break;
			}
			case 128: {
				MonoQueue6<128> mq;
				MONO<128>(mq);
				break;
			}
			default:
				throw std::runtime_error("missing monoqueue template");
		}
	}
}
BENCHMARK(BM_monoqueue_template)->RangeMultiplier(2)->Range(1<<1, 1<<7);












template <int k >
class MinSegmentTree {
	int i = 0;
	uint64_t segtree[2*k];

	constexpr int log2() {
		return std::ceil(std::log2(k));
	}
public:
	void set(uint32_t hash, uint32_t index) {
		int ind = i + k;
		if (++i == k) i = 0;

		segtree[ind] = (uint64_t)hash << 32 | ~index;
		for (int rep = 0; rep < log2(); rep++) {
			segtree[ind >> 1] = std::min(segtree[ind], segtree[ind ^ 1]);
			ind >>= 1;
		}
	}

	uint32_t min() {
		return ~segtree[1];
	}
};

template<uint32_t k>
void SEG(MinSegmentTree<k> &st) {
	uint32_t ind = 0;
	while(ind < k){
		st.set(hashes[ind], ind);
		ind++;
	}
	while(ind < INPUT_SIZE){
		st_outputs[ind-k] = st.min();
		st.set(hashes[ind], ind);
		ind++;
	}
	st_outputs[ind-k] = st.min();
}

static void BM_SegTree_template_unroll(benchmark::State& state){
    for(auto _ : state) {
		switch (state.range(0)) {
			case 2: {
				MinSegmentTree<2> st;
				SEG<2>(st);
				break;
			}
			case 4: {
				MinSegmentTree<4> st;
				SEG<4>(st);
				break;
			}
			case 8: {
				MinSegmentTree<8> st;
				SEG<8>(st);
				break;
			}
			case 16: {
				MinSegmentTree<16> st;
				SEG<16>(st);
				break;
			}
			case 32: {
				MinSegmentTree<32> st;
				SEG<32>(st);
				break;
			}
			case 64: {
				MinSegmentTree<64> st;
				SEG<64>(st);
				break;
			}
			case 128: {
				MinSegmentTree<128> st;
				SEG<128>(st);
				break;
			}
			default:
				throw std::runtime_error("undeclared k-mer seg tree");
		}

		benchmark::DoNotOptimize(st_outputs);
		benchmark::DoNotOptimize(st_outputs.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_SegTree_template_unroll)->RangeMultiplier(2)->Range(1<<1, 1<<7);

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


template<uint32_t k>
void SEG1() {
	ArraySegTree<k> st;
	int ind = 0;
	while(ind < k){
		st.set(hashes[ind], ind);
		ind++;
	}
	while(ind < INPUT_SIZE){
		st_outputs2[ind - k] = st.min();
		st.set(hashes[ind], ind);
		ind++;
	}
	st_outputs2[ind - k] = st.min();
}

static void BM_SegTree_array(benchmark::State& state){
    for(auto _ : state) {
		switch (state.range(0)) {
			case 2: {
				SEG1<2>();
				break;
			}
			case 4: {
				SEG1<4>();
				break;
			}
			case 8: {
				SEG1<8>();
				break;
			}
			case 16: {
				SEG1<16>();
				break;
			}
			case 32: {
				SEG1<32>();
				break;
			}
			case 64: {
				SEG1<64>();
				break;
			}
			case 128: {
				SEG1<128>();
				break;
			}
			default:
				throw std::runtime_error("undeclared k-mer array seg tree");
		}

		benchmark::DoNotOptimize(st_outputs2);
		benchmark::DoNotOptimize(st_outputs2.data());
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_SegTree_array)->RangeMultiplier(2)->Range(1<<1, 1<<7);

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
				st.set(hashes[ind], ind); \
				ind++; \
			} \
			while(ind < INPUT_SIZE){ \
				st_outputs5.push_back(st.min()); \
				st.set(hashes[ind], ind); \
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
// 			st.set(hashes[ind], ind);
// 			ind++;
// 		}
// 		while(ind < INPUT_SIZE){
// 			st_outputs2.push_back(st.min());
// 			st.set(hashes[ind], ind);
// 			ind++;
// 		}
// 		st_outputs2.push_back(st.min());
// 	}
// }
// BENCHMARK(BM_SegTree_128_unroll)->RangeMultiplier(2)->Range(1<<2, 1<<7);


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
//             st.set(hashes[ind], ind);
//             ind++;
//         }
//         while(ind < INPUT_SIZE){
//             st_outputs4.push_back(st.min());
//             st.set(hashes[ind], ind);
//             ind++;
//         }
// 		st_outputs4.push_back(st.min());
// 	}
// }
// BENCHMARK(BM_SegTree_128)->RangeMultiplier(2)->Range(1<<2, 1<<7);

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
//             st.set(hashes[ind], ind);
//             ind++;
//         }
//         while(ind < INPUT_SIZE){
//             st_outputs3.push_back(st.min());
//             st.set(hashes[ind], ind);
//             ind++;
//         }
// 		st_outputs3.push_back(st.min());
// 	}
// }
// BENCHMARK(BM_SegTree_pair)->RangeMultiplier(2)->Range(1<<2, 1<<7);



















static void BM_mset(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
		int k = state.range(0);
        int ind = 0;
        std::multiset<std::pair<uint32_t,int32_t>> mset;
        std::deque<decltype(mset)::iterator> vec;
        while(ind < k){
            vec.push_back(mset.emplace(hashes[ind], -ind));
            ind++;
        }
        int kick = 0;
        state.ResumeTiming();
        while(ind < INPUT_SIZE){
            mset_outputs[ind - k] = -mset.begin()->second;
            mset.erase(vec[kick]);
            vec[kick] = mset.emplace(hashes[ind], -ind);
            kick += 1;
            if(kick == k) kick = 0;
            ind++;
        }
		mset_outputs[ind - k] = -mset.begin()->second;
	}
}
BENCHMARK(BM_mset)->RangeMultiplier(2)->Range(1<<1, 1<<7);

static void BM_naive(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
		int k = state.range(0);
        state.ResumeTiming();

        for(int i = 0; i+k <= INPUT_SIZE; i++){
            uint64_t minAm = hashes[i];
			int ind = i;
            for(int j = i+1; j < i + k; j++){
				if (hashes[j] <= minAm) {
					minAm = hashes[j];
					ind = j;
				}
            }
            naive_outputs[i] = ind;
        }
	}
}
BENCHMARK(BM_naive)->RangeMultiplier(2)->Range(1<<1, 1<<7);


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
// 			if (hashes[i] < hashes[cur]) {
// 				cur = i;
// 			}
// 		}
//         for(int i = 0; i+k <= INPUT_SIZE; i++){
// 			if (hashes[i+k-1] < hashes[cur]) {
// 				cur = i+k-1;
// 			}
// 			if (cur == i-1) {
// 				cur++;
// 				for(int j = cur+1; j < i + k; j++){
// 					if (hashes[j] < hashes[cur]) {
// 						cur = j;
// 					}
// 				}
// 			}
//             naive_outputs2.push_back(cur);
//         }
// 	}
// }
// BENCHMARK(BM_naive2)->RangeMultiplier(2)->Range(1<<2, 1<<7);










int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    // sanity check
    assert(st_outputs == mset_outputs);
    assert(st_outputs == naive_outputs);
    assert(st_outputs == st_outputs2);
    assert(st_outputs == monoqueue_outputs6);

    std::cout << "Passed Asserts!" << std::endl;
}
