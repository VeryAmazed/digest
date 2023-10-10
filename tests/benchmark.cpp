//#include "nthash.hpp"
#include <digest/mod_minimizer.hpp>
#include <digest/window_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <fstream>
#include <benchmark/benchmark.h>

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

        digest::ModMin *dig = new digest::ModMin(str, 16, 17);
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
        digest::ModMin *dig = new digest::ModMin(str, 16, 17);
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

static void BM_WindowMinConstruction(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        state.ResumeTiming();

        digest::WindowMin *dig = new digest::WindowMin(str, 16, state.range(1));
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_WindowMinRoll(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        digest::WindowMin *dig = new digest::WindowMin(str, 16, state.range(1));
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

static void BM_SyncmerConstruction(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        state.ResumeTiming();

        digest::Syncmer *dig = new digest::Syncmer(str, 16, state.range(1));
        benchmark::DoNotOptimize(dig);
        
        state.PauseTiming();
        delete dig;
        state.ResumeTiming();
    }
}


static void BM_SyncmerRoll(benchmark::State& state){
    for(auto _ : state){
        state.PauseTiming();
        std::string str = bench_strs[0].substr(0, state.range(0));
        digest::Syncmer *dig = new digest::Syncmer(str, 16, state.range(1));
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


BENCHMARK(BM_ModMinConstruction)->RangeMultiplier(10)->Range(1e3, 1e6)->Complexity();
BENCHMARK(BM_ModMinRoll)->RangeMultiplier(10)->Range(1e3, 1e6)->Complexity();

/*
BENCHMARK(BM_WindowMinConstruction)->ArgsProduct({
    benchmark::CreateRange(1e3, 1e6, 10),
      benchmark::CreateRange(8, 64, 2)});
BENCHMARK(BM_WindowMinRoll)->ArgsProduct({
    benchmark::CreateRange(1e3, 1e6, 10),
      benchmark::CreateRange(8, 64, 2)});

BENCHMARK(BM_SyncmerConstruction)->ArgsProduct({
    benchmark::CreateRange(1e3, 1e6, 10),
      benchmark::CreateRange(8, 64, 2)});
BENCHMARK(BM_SyncmerRoll)->ArgsProduct({
    benchmark::CreateRange(1e3, 1e6, 10),
      benchmark::CreateRange(8, 64, 2)});
*/

int main(int argc, char** argv)
{
   setupStrings();
   ::benchmark::Initialize(&argc, argv);
   ::benchmark::RunSpecifiedBenchmarks();
}