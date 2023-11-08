// perf record --call-graph dwarf bench
// perf report -g

#include <digest/mod_minimizer.hpp>
#include <digest/window_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <fstream>
#include <benchmark/benchmark.h>
#include <nthash/nthash.hpp>

#define DEFAULT_LARGE_WIND 16
#define DEFAULT_KMER_LEN 16
#define DEFAULT_STR_LEN 1e4

std::vector<std::string> bench_strs;
std::string s;
std::string s1;
std::string s2;


void setupStrings(){
	std::string files[] = {
		"../tests/benchmark_strings/ACTG.txt",
	};
	
	for (auto& file : files) {
		std::ifstream ifs(file);
		ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        std::string str;
		ifs >> str;
		bench_strs.push_back(str);
	}
}

static void random(const benchmark::State& state) {
	s = bench_strs[0].substr(0, state.range(0));
}

static void ecoli(const benchmark::State& state) {
	s = bench_strs[1].substr(0, state.range(0));
}

static void random_append_seq(const benchmark::State& state){
	s1 = bench_strs[0].substr(0, state.range(0));
	s2 = s = bench_strs[0].substr(state.range(0), state.range(1));
}

// static void BM_NtHashConstruction(benchmark::State& state){
// 	for(auto _ : state) {
// 		nthash::NtHash dig(s, 1, DEFAULT_KMER_LEN);
// 		benchmark::DoNotOptimize(dig);
// 		benchmark::ClobberMemory();
// 	}
// 	state.SetComplexityN(state.range(0));
// }
// BENCHMARK(BM_NtHashConstruction)->Range(1<<6, 1<<18)->Setup(random)->Complexity();

static void BM_NtHashRoll(benchmark::State& state) {
	for(auto _ : state) {
		nthash::NtHash dig(s, 1, DEFAULT_KMER_LEN);

		while (dig.roll())
			benchmark::DoNotOptimize(*dig.hashes());
	}
	state.SetComplexityN(state.range(0));
}
// BENCHMARK(BM_NtHashRoll)->Range(1<<6, 1<<18)->Setup(random)->Complexity();


// static void BM_ModMinConstruction(benchmark::State& state){
// 	for(auto _ : state) {
// 		digest::ModMin dig(s, DEFAULT_KMER_LEN, 17);
// 		benchmark::DoNotOptimize(dig);
// 		benchmark::ClobberMemory();
// 	}
// 	state.SetComplexityN(state.range(0));
// }
// BENCHMARK(BM_ModMinConstruction)->Range(1<<6, 1<<18)->Setup(random)->Complexity();

static void BM_ModMinRoll(benchmark::State& state) {
	for(auto _ : state) {
		digest::ModMin dig(s, DEFAULT_KMER_LEN, 17);
		std::vector<size_t> vec;
		vec.reserve(state.range(0));

		benchmark::DoNotOptimize(vec);
		dig.roll_minimizer(state.range(0), vec);
		benchmark::ClobberMemory();
	}
	state.SetComplexityN(state.range(0));
}
// BENCHMARK(BM_ModMinRoll)->Range(1<<6, 1<<18)->Setup(random)->Complexity();


// static void BM_WindowMinConstructionFixWind(benchmark::State& state){
//     for(auto _ : state){
//         digest::WindowMin dig(s, DEFAULT_KMER_LEN, DEFAULT_LARGE_WIND);
//         benchmark::DoNotOptimize(dig);
// 		benchmark::ClobberMemory();
//     }
// 	state.SetComplexityN(state.range(0));
// }
// BENCHMARK(BM_WindowMinConstructionFixWind)->Range(1<<6, 1<<18)->Setup(random)->Complexity();

static void BM_WindowMinRollFixWind(benchmark::State& state){
    for(auto _ : state){
        digest::WindowMin dig(s, DEFAULT_KMER_LEN,DEFAULT_LARGE_WIND);
        std::vector<size_t> vec;
        vec.reserve(state.range(0));

		benchmark::DoNotOptimize(vec);
        dig.roll_minimizer(state.range(0), vec);
		benchmark::ClobberMemory();
    }
	state.SetComplexityN(state.range(0));
}
// BENCHMARK(BM_WindowMinRollFixWind)->Range(1<<6, 1<<18)->Setup(random)->Complexity();




// static void BM_WindowMinConstructionFixLen(benchmark::State& state){
//     for(auto _ : state){
//         digest::WindowMin dig(s, DEFAULT_KMER_LEN, state.range(0));
//         benchmark::DoNotOptimize(dig);
// 		benchmark::ClobberMemory();
//     }
// 	state.SetComplexityN(state.range(0));
// }
// BENCHMARK(BM_WindowMinConstructionFixLen)->Range(1<<6, 1<<18)->Setup(random)->Complexity();

static void BM_WindowMinRollFixLen(benchmark::State& state){
    for(auto _ : state){
        digest::WindowMin dig(s, DEFAULT_KMER_LEN, state.range(0));
        std::vector<size_t> vec;
        vec.reserve(DEFAULT_STR_LEN);

		benchmark::DoNotOptimize(vec);
        dig.roll_minimizer(DEFAULT_STR_LEN, vec);
		benchmark::ClobberMemory();
    }
	state.SetComplexityN(state.range(0));
}
// BENCHMARK(BM_WindowMinRollFixLen)->Range(1<<6, 1<<18)->Setup(random)->Complexity();




// static void BM_SyncmerConstructionFixWind(benchmark::State& state){
//     for(auto _ : state){
//         digest::Syncmer dig(s, DEFAULT_KMER_LEN, DEFAULT_LARGE_WIND);
//         benchmark::DoNotOptimize(dig);
// 		benchmark::ClobberMemory();
//     }
// 	state.SetComplexityN(state.range(0));
// }
// BENCHMARK(BM_SyncmerConstructionFixWind)->Range(1<<6, 1<<18)->Setup(random)->Complexity();


static void BM_SyncmerRollFixWind(benchmark::State& state){
    for(auto _ : state){
        digest::Syncmer dig(s, DEFAULT_KMER_LEN,DEFAULT_LARGE_WIND);
        std::vector<size_t> vec;
        vec.reserve(state.range(0));

		benchmark::DoNotOptimize(vec);
        dig.roll_minimizer(state.range(0), vec);
		benchmark::ClobberMemory();
    }
	state.SetComplexityN(state.range(0));
}
// BENCHMARK(BM_SyncmerRollFixWind)->Range(1<<6, 1<<18)->Setup(random)->Complexity();


// static void BM_SyncmerConstructionFixLen(benchmark::State& state){
//     for(auto _ : state){
//         digest::Syncmer dig(s, DEFAULT_KMER_LEN, state.range(0));
//         benchmark::DoNotOptimize(dig);
// 		benchmark::ClobberMemory();
//     }
// 	state.SetComplexityN(state.range(0));
// }
// BENCHMARK(BM_SyncmerConstructionFixLen)->Range(1<<6, 1<<18)->Setup(random)->Complexity();


static void BM_SyncmerRollFixLen(benchmark::State& state){
    for(auto _ : state){
        digest::Syncmer dig(s, DEFAULT_KMER_LEN, state.range(0));
        std::vector<size_t> vec;
        vec.reserve(DEFAULT_STR_LEN);

		benchmark::DoNotOptimize(vec);
        dig.roll_minimizer(DEFAULT_STR_LEN, vec);
		benchmark::ClobberMemory();
    }
	state.SetComplexityN(state.range(0));
}
// BENCHMARK(BM_SyncmerRollFixLen)->Range(1<<6, 1<<18)->Setup(random)->Complexity();

// more of a sanity check
static void BM_append_seq(benchmark::State& state){
	int kmer_len = 64;
	for(auto _ : state) {
		digest::ModMin dig(s1, kmer_len, 17);
		dig.append_seq(s2);
		std::vector<size_t> vec;
		vec.reserve(state.range(0) + state.range(1));

		benchmark::DoNotOptimize(vec);
		dig.roll_minimizer(state.range(0), vec);
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_ModMinRoll)->Arg(127)->Setup(random);
BENCHMARK(BM_append_seq)->Args({63, 64})->Setup(random_append_seq);

int main(int argc, char** argv)
{
   setupStrings();
   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
}
