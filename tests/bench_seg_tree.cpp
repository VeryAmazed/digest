#include <benchmark/benchmark.h>
#include <cwchar>
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
std::vector<int> monoqueue_outputs;
std::vector<int> monoqueue_outputs2;

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
        state.ResumeTiming();

        int ind = 0;
        MinSegmentTree<int> st(k);
        while(ind < k){
            st.set(ind, inputs[ind]);
            ind++;
        }
        int kick = 0;
        while(ind < INPUT_SIZE){
            st_outputs.push_back(st.segtree[1]);
            st.set(kick, inputs[ind]);
            
            kick += 1;
            if(kick == k) kick = 0;
            
            //kick = (kick+1)%k;
            ind++;
        }
		st_outputs.push_back(st.segtree[1]);
	}
}
BENCHMARK(BM_SegTree)->RangeMultiplier(2)->Range(1<<2, 1<<6);

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
		mset_outputs.push_back(*(mset.begin()));
	}
}
//BENCHMARK(BM_mset)->RangeMultiplier(2)->Range(1<<2, 1<<6);

static void BM_naive(benchmark::State& state){
    for(auto _ : state) {
        state.PauseTiming();
        naive_outputs.clear();
        naive_outputs.reserve(INPUT_SIZE);
		int k = state.range(0);

        state.ResumeTiming();
        for(int i = 0; i+k <= INPUT_SIZE; i++){
            int minAm = inputs[i];
            for(int j = i+1; j < i + k; j++){
                minAm = std::min(minAm, inputs[j]);
            }
            naive_outputs.push_back(minAm);
        }
	}
}
BENCHMARK(BM_naive)->RangeMultiplier(2)->Range(1<<2, 1<<6);









class MonoQueue {
	std::pair<int,int> deq[64];
	int head = 0, tail = 0;
	const int k;

	// cannot call this while full
	bool empty() {
		return head == tail;
	}

	void emplace_back(int l, int r) {
		deq[tail] = {l, r};
		if (tail == k-1) tail = 0;
		else tail++;
	}

	void pop_front() {
		if (head == k-1) head = 0;
		else head++;
	}

	void pop_back() {
		if (--tail == -1) tail = k - 1;
	}

	int front1() {
		return deq[head].first;
	}

	int back() {
		return deq[tail == 0 ? k - 1 : tail - 1].second;
	}
public:

	MonoQueue(int k) : k(k) {}

	void add (int i) {
		while (!empty() && back() >= inputs[i]) {
			pop_back();
		}
		if (!empty() && front1() == i - k) {
			pop_front();
		}
		emplace_back(i, inputs[i]);
	}

	int min() {
		return deq[head].second;
	}
};







class MonoQueue2 {
	std::pair<int,int> deq[64];
	int head = 0, tail = 0, k;

	bool empty() {
		return head == tail;
	}

	void emplace_back(int l, int r) {
		deq[tail] = {l, r};
		if (++tail == k) tail = 0;
	}

	void pop_front() {
		if (++head == k) head = 0;
	}

	void pop_back() {
		if (--tail == -1) tail = k - 1;
	}

	int front1() {
		return deq[head].first;
	}

	int back() {
		return deq[tail == 0 ? k - 1 : tail - 1].second;
	}

public:
	MonoQueue2(int k) : k(k) {}

	void add (int i) {
		while (!empty() and back() >= inputs[i]) {
			pop_back();
		}
		if (!empty() and front1() == i - k) {
			pop_front();
		}
		emplace_back(i, inputs[i]);
	}

	int min() {
		return deq[head].second;
	}
};







static void BM_monoqueue(benchmark::State& state) {
    for(auto _ : state) {
        state.PauseTiming();
        monoqueue_outputs.clear();
		int k = state.range(0);
        state.ResumeTiming();

		MonoQueue mq(k);

		int i = 0;
		while (i < k) {
			mq.add(i++);
		}

        for (; i < INPUT_SIZE; i++){
			monoqueue_outputs.push_back(mq.min());
			mq.add(i);
        }
		monoqueue_outputs.push_back(mq.min());
	}
}
BENCHMARK(BM_monoqueue)->RangeMultiplier(2)->Range(1<<2, 1<<6);

static void BM_monoqueue2(benchmark::State& state) {
    for(auto _ : state) {
        state.PauseTiming();
        monoqueue_outputs2.clear();
		int k = state.range(0);
        state.ResumeTiming();

		MonoQueue2 mq(k);

		int i = 0;
		while (i < k) {
			mq.add(i++);
		}

        for (; i < INPUT_SIZE; i++){
			monoqueue_outputs2.push_back(mq.min());
			mq.add(i);
        }
		monoqueue_outputs2.push_back(mq.min());
	}
}
BENCHMARK(BM_monoqueue2)->RangeMultiplier(2)->Range(1<<2, 1<<6);





int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    // sanity check
    // if it's correct for 1 window size, it's probably correct for all
    //assert(st_outputs == mset_outputs);
    assert(st_outputs == naive_outputs);
    assert(st_outputs == monoqueue_outputs);
    assert(st_outputs == monoqueue_outputs2);
    std::cout << "Passed Asserts!" << std::endl;
}
