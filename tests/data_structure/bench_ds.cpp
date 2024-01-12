#include <digest/data_structure.hpp>
#include <benchmark/benchmark.h>

#include <cstdint>
#include <iostream>
#include <random>

// segtee wins at 12
// naive2 wins at 17

const int INPUT_SIZE = 1e7;
std::array<uint32_t,2*INPUT_SIZE> hashes;

std::map<int,std::map<int,std::array<uint32_t, INPUT_SIZE>>> all;

void setupInput(){
	std::random_device rd; // seed
	std::mt19937 gen(rd()); // generator
	std::uniform_int_distribution<uint32_t> distrib(0,UINT32_MAX); // [0, 2**32]
	for (uint32_t &h : hashes) {
        h = distrib(gen);
    }

	// edge test for ties
	// for (int i = 0; i < 2*INPUT_SIZE; i++) {
	// 	hashes[i] = hashes[0];
	// }
}


template<int k, class T, int out>
static void BM(benchmark::State& state){
	auto &temp = all[out][k];
    for(auto _ : state) {
		T ds(k);
		for (int i = 0; i < k-1; i++) {
			ds.insert(i, hashes[i]);
		}
        for (int i = 0; i < INPUT_SIZE; i++) {
            temp[i] = ds.insert(i+k-1, hashes[i+k-1]);
        }
		benchmark::ClobberMemory();
	}
}

#define test(name, out) \
	BENCHMARK_TEMPLATE(BM, 2, name<2>, out); \
	BENCHMARK_TEMPLATE(BM, 3, name<3>, out); \
	BENCHMARK_TEMPLATE(BM, 4, name<4>, out); \
	// BENCHMARK_TEMPLATE(BM, 5, name<5>, out); \
	// BENCHMARK_TEMPLATE(BM, 8, name<8>, out); \
	// BENCHMARK_TEMPLATE(BM, 9, name<9>, out); \
	// BENCHMARK_TEMPLATE(BM, 12, name<12>, out); \
	// BENCHMARK_TEMPLATE(BM, 16, name<16>, out); \
	// BENCHMARK_TEMPLATE(BM, 17, name<17>, out); \
	// BENCHMARK_TEMPLATE(BM, 32, name<32>, out); \
	// BENCHMARK_TEMPLATE(BM, 33, name<33>, out); \
	// BENCHMARK_TEMPLATE(BM, 64, name<64>, out); \
	// BENCHMARK_TEMPLATE(BM, 96, name<96>, out); \
	// BENCHMARK_TEMPLATE(BM, 128, name<128>, out); \
	// BENCHMARK_TEMPLATE(BM, 256, name<256>, out); \
	// BENCHMARK_TEMPLATE(BM, 512, name<512>, out); \
	// BENCHMARK_TEMPLATE(BM, 1024, name<1024>, out); \

#define test2(name, out) \
	BENCHMARK_TEMPLATE(BM, 2, name, out); \
	BENCHMARK_TEMPLATE(BM, 3, name, out); \
	BENCHMARK_TEMPLATE(BM, 4, name, out); \
	// BENCHMARK_TEMPLATE(BM, 5, name, out); \
	// BENCHMARK_TEMPLATE(BM, 8, name, out); \
	// BENCHMARK_TEMPLATE(BM, 9, name, out); \
	// BENCHMARK_TEMPLATE(BM, 12, name, out); \
	// BENCHMARK_TEMPLATE(BM, 16, name, out); \
	// BENCHMARK_TEMPLATE(BM, 17, name, out); \
	// BENCHMARK_TEMPLATE(BM, 32, name, out); \
	// BENCHMARK_TEMPLATE(BM, 33, name, out); \
	// BENCHMARK_TEMPLATE(BM, 64, name, out); \
	// BENCHMARK_TEMPLATE(BM, 96, name, out); \
	// BENCHMARK_TEMPLATE(BM, 128, name, out); \
	// BENCHMARK_TEMPLATE(BM, 256, name, out); \
	// BENCHMARK_TEMPLATE(BM, 512, name, out); \
	// BENCHMARK_TEMPLATE(BM, 1024, name, out); \

test(Naive, 0);
test(Naive2, 1);
test(MonoQueue, 2);
test(SegmentTree, 3);
test(Set, 4);
test2(Adaptive, 5);

int main(int argc, char** argv)
{
    setupInput();
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    // sanity check
	for (auto &[_, m] : all) {
		assert(m == all[0]);
	}

    std::cout << "Passed Asserts!" << std::endl;
}
