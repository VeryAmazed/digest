//#include "nthash.hpp"
#include <digest/mod_minimizer.hpp>
#include <digest/window_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <fstream>
#include <benchmark/benchmark.h>

#define DEFAULT_LARGE_WIND 16
#define DEFAULT_KMER_LEN 16
#define DEFAULT_STR_LEN 1e4

std::vector<std::string> bench_strs;

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

static void BM_ModMinConstruction(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        state.ResumeTiming();

        digest::ModMin *dig = new digest::ModMin(str, DEFAULT_KMER_LEN, 17);
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_ModMinRoll(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        digest::ModMin *dig = new digest::ModMin(str, DEFAULT_KMER_LEN, 17);
        std::vector<size_t> *vec = new std::vector<size_t>();
        vec->reserve(state.range(0));
        benchmark::DoNotOptimize(vec);
        state.ResumeTiming();
        dig->roll_minimizer(state.range(0), *vec);
        state.PauseTiming();
        delete vec;
        delete dig;
        state.ResumeTiming();
    }
}

static void BM_WindowMinConstructionFixWind(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        state.ResumeTiming();

        digest::WindowMin *dig = new digest::WindowMin(str, DEFAULT_KMER_LEN, DEFAULT_LARGE_WIND);
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_WindowMinRollFixWind(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        digest::WindowMin *dig = new digest::WindowMin(str, DEFAULT_KMER_LEN,DEFAULT_LARGE_WIND);
        std::vector<size_t> *vec = new std::vector<size_t>();
        vec->reserve(state.range(0));
        benchmark::DoNotOptimize(vec);
        state.ResumeTiming();
        dig->roll_minimizer(state.range(0), *vec);
        state.PauseTiming();
        delete vec;
        delete dig;
        state.ResumeTiming();
    }
}

static void BM_WindowMinConstructionFixLen(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, DEFAULT_STR_LEN);
        state.ResumeTiming();

        digest::WindowMin *dig = new digest::WindowMin(str, DEFAULT_KMER_LEN, state.range(0));
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_WindowMinRollFixLen(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, DEFAULT_STR_LEN);
        digest::WindowMin *dig = new digest::WindowMin(str, DEFAULT_KMER_LEN, state.range(0));
        std::vector<size_t> *vec = new std::vector<size_t>();
        vec->reserve(DEFAULT_STR_LEN);
        benchmark::DoNotOptimize(vec);
        state.ResumeTiming();
        dig->roll_minimizer(DEFAULT_STR_LEN, *vec);
        state.PauseTiming();
        delete vec;
        delete dig;
        state.ResumeTiming();
    }
}

static void BM_SyncmerConstructionFixWind(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        state.ResumeTiming();

        digest::Syncmer *dig = new digest::Syncmer(str, DEFAULT_KMER_LEN, DEFAULT_LARGE_WIND);
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_SyncmerRollFixWind(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        digest::Syncmer *dig = new digest::Syncmer(str, DEFAULT_KMER_LEN,DEFAULT_LARGE_WIND);
        std::vector<size_t> *vec = new std::vector<size_t>();
        vec->reserve(state.range(0));
        benchmark::DoNotOptimize(vec);
        state.ResumeTiming();
        dig->roll_minimizer(state.range(0), *vec);
        state.PauseTiming();
        delete vec;
        delete dig;
        state.ResumeTiming();
    }
}

static void BM_SyncmerConstructionFixLen(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, DEFAULT_STR_LEN);
        state.ResumeTiming();

        digest::Syncmer *dig = new digest::Syncmer(str, DEFAULT_KMER_LEN, state.range(0));
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_SyncmerRollFixLen(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, DEFAULT_STR_LEN);
        digest::Syncmer *dig = new digest::Syncmer(str, DEFAULT_KMER_LEN, state.range(0));
        std::vector<size_t> *vec = new std::vector<size_t>();
        vec->reserve(DEFAULT_STR_LEN);
        benchmark::DoNotOptimize(vec);
        state.ResumeTiming();
        dig->roll_minimizer(DEFAULT_STR_LEN, *vec);
        state.PauseTiming();
        delete vec;
        delete dig;
        state.ResumeTiming();
    }
}

BENCHMARK(BM_ModMinConstruction)->RangeMultiplier(10)->Range(1e3, 1e6);
BENCHMARK(BM_ModMinRoll)->RangeMultiplier(10)->Range(1e3, 1e6);

BENCHMARK(BM_WindowMinConstructionFixWind)->RangeMultiplier(10)->Range(1e3, 1e6);
BENCHMARK(BM_WindowMinRollFixWind)->RangeMultiplier(10)->Range(1e3, 1e6);

BENCHMARK(BM_WindowMinConstructionFixLen)->RangeMultiplier(2)->Range(8, 64);
BENCHMARK(BM_WindowMinRollFixLen)->RangeMultiplier(2)->Range(8, 64);

BENCHMARK(BM_SyncmerConstructionFixWind)->RangeMultiplier(10)->Range(1e3, 1e6);
BENCHMARK(BM_SyncmerRollFixWind)->RangeMultiplier(10)->Range(1e3, 1e6);

BENCHMARK(BM_SyncmerConstructionFixLen)->RangeMultiplier(2)->Range(8, 64);
BENCHMARK(BM_SyncmerRollFixLen)->RangeMultiplier(2)->Range(8, 64);
int main(int argc, char** argv)
{
   setupStrings();
   ::benchmark::Initialize(&argc, argv);
   ::benchmark::RunSpecifiedBenchmarks();
}