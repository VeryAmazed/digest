#include <benchmark/benchmark.h>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <deque>

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
		for (; ind > 1; ind /= 2) {
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
		int k = state.range(0);
        int ind = 0;
        MinSegmentTree<int> st(k);
        while(ind < k){
            st.set(ind, inputs[ind]);
            ind++;
        }
        int kick = 0;
        state.ResumeTiming();
        int a = 0;
        while(ind < INPUT_SIZE){
            benchmark::DoNotOptimize(a);
            a = st.segtree[1];
            benchmark::ClobberMemory();
            st.set(kick, inputs[ind]);
            kick = (kick+1)%k;
            ind++;
        }
	}
}
BENCHMARK(BM_SegTree)->RangeMultiplier(2)->Range(1<<2, 1<<6);

static void BM_mset_dq(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
		int k = state.range(0);
        int ind = 0;
        std::multiset<int> mset;
        std::deque<std::multiset<int>::iterator> dq;
        benchmark::DoNotOptimize(mset);
        benchmark::DoNotOptimize(dq);
        while(ind < k){
            dq.push_back(mset.insert(inputs[ind]));
            ind++;
        }
        state.ResumeTiming();

        while(ind < INPUT_SIZE){
            mset.erase(dq.front());
            dq.pop_front();
            dq.push_back(mset.insert(inputs[ind]));

            ind++;
        }
	}
}
BENCHMARK(BM_mset_dq)->RangeMultiplier(2)->Range(1<<2, 1<<6);

int main(int argc, char** argv)
{
   setupInput();
   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
}
