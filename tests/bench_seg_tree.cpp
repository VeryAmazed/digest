#include <benchmark/benchmark.h>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <deque>
#include <cassert>
#include <iostream>
#include <algorithm>

#define INPUT_SIZE 1e4

template <class T> class MinSegmentTree {
  public:
	/** The operation to use for combining two elements. (Must be associative)
	 */
	T comb(T a, T b) { return std::min(a, b); }
	const T DEFAULT = 2e9;  // Will overflow if T is an int

	std::vector<T> segtree;
	int len;

	MinSegmentTree(int len) : segtree(len * 2, DEFAULT), len(len) {}

	/** Sets the value at ind to val. */
	void set(int ind, T val) {
		//assert(0 <= ind && ind < len);
		ind += len;
		segtree[ind] = val;
		for (; ind > 1; ind >>= 1) {
			segtree[ind >> 1] = comb(segtree[ind], segtree[ind ^ 1]);
		}
	}

	/** @return the minimum element in the range [start, end) */
	T range_min(int start, int end) {
		//assert(0 <= start && start < len && 0 < end && end <= len);
		T sum = DEFAULT;
		for (start += len, end += len; start < end; start /= 2, end /= 2) {
			if ((start & 1) != 0) { sum = comb(sum, segtree[start++]); }
			if ((end & 1) != 0) { sum = comb(sum, segtree[--end]); }
		}
		return sum;
	}
};

std::vector<int> inputs;
std::vector<int> st_outputs;
std::vector<int> mset_outputs;
std::vector<int> naive_outputs;

void setupInput(){
    std::string path = "../tests/benchmark_seg_tree_input/nums.txt";
    std::ifstream ifs(path);
    ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    for(int i =0; i < INPUT_SIZE; i++){
        int a;
        ifs >> a;
        inputs.push_back(a);
    }
}

static void BM_SegTree(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
        st_outputs.clear();
        st_outputs.reserve(INPUT_SIZE);
		int k = state.range(0);
        int ind = 0;
        MinSegmentTree<int> st(k);
        while(ind < k){
            st.set(ind, inputs[ind]);
            ind++;
        }
        int kick = 0;
        state.ResumeTiming();
        while(ind < INPUT_SIZE){
            st_outputs.push_back(st.segtree[1]);
            st.set(kick, inputs[ind]);
            kick += 1;
            if(kick == k) kick = 0;
            ind++;
        }
	}
}
BENCHMARK(BM_SegTree)->RangeMultiplier(2)->Range(1<<2, 1<<8);

static void BM_mset(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
        mset_outputs.clear();
        mset_outputs.reserve(INPUT_SIZE);
		int k = state.range(0);
        int ind = 0;
        std::multiset<int> mset;
        std::deque<std::multiset<int>::iterator> vec;
        while(ind < k){
            vec.push_back(mset.insert(inputs[ind]));
            ind++;
        }
        int kick = 0;
        state.ResumeTiming();
        while(ind < INPUT_SIZE){
            mset_outputs.push_back(*(mset.begin()));
            mset.erase(vec[kick]);
            vec[kick] = mset.insert(inputs[ind]);
            kick += 1;
            if(kick == k) kick = 0;
            ind++;
        }
	}
}
BENCHMARK(BM_mset)->RangeMultiplier(2)->Range(1<<2, 1<<8);

static void BM_naive(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
        naive_outputs.clear();
        naive_outputs.reserve(INPUT_SIZE);
		int k = state.range(0);

        state.ResumeTiming();
        for(int i = 0; i+k < INPUT_SIZE; i++){
            int minAm = inputs[i];
            for(int j = i+1; j < i + k; j++){
                minAm = std::min(minAm, inputs[j]);
            }
            naive_outputs.push_back(minAm);
        }
	}
}
BENCHMARK(BM_naive)->RangeMultiplier(2)->Range(1<<2, 1<<8);

int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    // sanity check
    // if it's correct for 1 window size, it's probably correct for all
    assert(st_outputs == mset_outputs);
    assert(st_outputs == naive_outputs);
    std::cout << "Passed Asserts!" << std::endl;
}
