#include <catch2/catch_test_macros.hpp>
#include "digest/thread_out.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

std::vector<std::string> test_strs;

void setupStrings(){
	std::string files[] = {
		"../tests/test/A.txt",
		"../tests/test/a_lowercase.txt",
		"../tests/test/salmonella_enterica.txt",
		"../tests/test/salmonella_lowercase.txt",
		"../tests/test/random.txt",
		"../tests/test/random_lowercase.txt",
		"../tests/test/N.txt",
	};
	
	for (auto& file : files) {
		std::ifstream ifs(file);
		ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		std::string str;
		ifs >> str;
		test_strs.push_back(str);
	}
}

template <class T>
std::vector<T> multi_to_single_vec(std::vector<std::vector<T>> vec){
    std::vector<T> ret_vec;
    for(size_t i = 0; i < vec.size(); i++){
        for(auto a : vec[i]){
            ret_vec.push_back(a);
        }
    }
    return ret_vec;
}

void test_thread_mod(unsigned thread_count,
    std::string str, unsigned k, uint64_t mod, uint64_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
    std::vector<uint32_t> single_thread;
	std::vector<std::vector<uint32_t>> vec;
    digest::ModMin dig(str, k, mod, congruence, start, minimized_h);
    dig.roll_minimizer(str.size(), single_thread);
    thread_out::thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h);
    std::vector<uint32_t> multi_thread = multi_to_single_vec(vec);
    
    REQUIRE(single_thread.size() == multi_thread.size());
	for(size_t i = 0; i < single_thread.size(); i++){
		CHECK(single_thread[i] == multi_thread[i]);
	}
}

void test_thread_wind(unsigned thread_count,
    std::string str, unsigned k, unsigned large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
    std::vector<uint32_t> single_thread;
	std::vector<std::vector<uint32_t>> vec;
    digest::WindowMin<data_structure::Adaptive> dig(str, k, large_wind_kmer_am, start, minimized_h);
    dig.roll_minimizer(str.size(), single_thread);
    thread_out::thread_wind<data_structure::Adaptive>(thread_count, vec, str, k, large_wind_kmer_am, start, minimized_h);
    std::vector<uint32_t> multi_thread = multi_to_single_vec(vec);
    
    REQUIRE(single_thread.size() == multi_thread.size());
	for(size_t i = 0; i < single_thread.size(); i++){
		CHECK(single_thread[i] == multi_thread[i]);
	}
}

void test_thread_sync(unsigned thread_count,
    std::string str, unsigned k, unsigned large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
    std::vector<uint32_t> single_thread;
	std::vector<std::vector<uint32_t>> vec;
    digest::Syncmer<data_structure::Adaptive> dig(str, k, large_wind_kmer_am, start, minimized_h);
    dig.roll_minimizer(str.size(), single_thread);
    thread_out::thread_sync<data_structure::Adaptive>(thread_count, vec, str, k, large_wind_kmer_am, start, minimized_h);
    std::vector<uint32_t> multi_thread = multi_to_single_vec(vec);
    
    REQUIRE(single_thread.size() == multi_thread.size());
	for(size_t i = 0; i < single_thread.size(); i++){
		CHECK(single_thread[i] == multi_thread[i]);
	}
}

TEST_CASE("thread_mod function testing"){
    setupStrings();
    SECTION("Throw Errors"){
        std::string str = "ACTGACTG";
        unsigned thread_count = 4;
        std::vector<std::vector<uint32_t>> vec;
        unsigned k = 4;
        uint64_t mod = 1e9+7;
        uint64_t congruence = 0;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // k < 4
        k = 3;
        CHECK_THROWS_AS(thread_out::thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h), thread_out::BadThreadOutParams);
        k = 4;
        // start >= len
        start = str.size();
        CHECK_THROWS_AS(thread_out::thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h), thread_out::BadThreadOutParams);
        start = 0;
        // num_kmers is negative
        start = 7;
        CHECK_THROWS_AS(thread_out::thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h), thread_out::BadThreadOutParams);
        start = 0;
        // num_kmers < thread_count
        thread_count = 6;
        CHECK_THROWS_AS(thread_out::thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h), thread_out::BadThreadOutParams);
        thread_count = 4;
    }

    SECTION("Special Cases"){
        unsigned thread_count = 4;
        unsigned k = 4;
        uint64_t mod = 3;
        uint64_t congruence = 0;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // only 1 thread
        thread_count = 1;
        for(int i = 0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_mod(thread_count, str, k, mod, congruence, start, minimized_h);
            }
        }

        // each thread gets 1 kmer
        thread_count = 96;
        for(int i =0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_mod(thread_count, str, k, mod, congruence, start, minimized_h);
            }
        }

        // some threads get 2 kmers, the rest get 1
        thread_count = 50;
        for(int i = 0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_mod(thread_count, str, k, mod, congruence, start, minimized_h);
            }
        }
    }

    SECTION("Full Testing"){
        unsigned thread_count = 4;
        unsigned k = 4;
        uint64_t mod = 17;
        uint64_t congruence = 0;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // the string changes
        // thread_count changes
        // start changes
        for(int i =0; i < 4; i++){
            for(int i = 0; i < 4; i += 2){
                for(int j = 4; j <= 64; j += 4){
                    thread_count = j;
                    for(int l = 0; l <= 96; l += 13){
                        start = l;
                        test_thread_mod(thread_count, test_strs[i], k, mod, congruence, start, minimized_h);
                    }
                }
            }
        }
    }
}

TEST_CASE("thread_wind function testing"){
    setupStrings();
    SECTION("Throw Errors"){
        std::string str = "ACTGACTGACTG";
        unsigned thread_count = 4;
        std::vector<std::vector<uint32_t>> vec;
        unsigned k = 4;
        const uint32_t large_wind_kmer_am = 4;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // large_wind_kmer_am is 0
        // CHECK_THROWS_AS(thread_out::thread_wind<data_>(thread_count, vec, str, k, start, minimized_h), thread_out::BadThreadOutParams);
        // num_lwinds is negative
        start = 9;
        CHECK_THROWS_AS(thread_out::thread_wind<data_structure::Adaptive>(thread_count, vec, str, k, large_wind_kmer_am, start, minimized_h), thread_out::BadThreadOutParams);
        start = 0;
        // num_lwinds < thread_count
        thread_count = 8;
        CHECK_THROWS_AS(thread_out::thread_wind<data_structure::Adaptive>(thread_count, vec, str, k, large_wind_kmer_am, start, minimized_h), thread_out::BadThreadOutParams);
        thread_count = 4;
    }

    SECTION("Special Cases"){
        unsigned thread_count = 4;
        unsigned k = 4;
        const uint32_t large_wind_kmer_am = 8; 
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        
        // only 1 thread
        thread_count = 1;
        for(int i = 0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_wind(thread_count, str, k, large_wind_kmer_am, start, minimized_h);
            }
        }
        
        // each thread gets 1 lwind
        thread_count = 86;
        for(int i =0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_wind(thread_count, str, k, large_wind_kmer_am, start, minimized_h);
            }
        }
        
        // some threads get 2 kmers, the rest get 1
        thread_count = 50;
        for(int i = 0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_wind(thread_count, str, k, large_wind_kmer_am, start, minimized_h);
            }
        }
        
    }
    
    SECTION("Full Testing"){
        unsigned thread_count = 4;
        unsigned k = 4;
        const uint32_t large_wind_kmer_am = 8;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // the string changes
        // thread_count changes
        // start changes
        for(int i =0; i < 4; i++){
            for(int i = 0; i < 4; i += 2){
                for(int j = 4; j <= 64; j += 4){
                    thread_count = j;
                    for(int l = 0; l <= 96; l += 13){
                        start = l;
                        test_thread_wind(thread_count, test_strs[i], k, large_wind_kmer_am, start, minimized_h);
                    }
                }
            }
        }
    }
}

TEST_CASE("thread_sync function testing"){
    setupStrings();
    SECTION("Special Cases"){
        unsigned thread_count = 4;
        unsigned k = 4;
        const uint32_t large_wind_kmer_am = 8; 
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        
        // only 1 thread
        thread_count = 1;
        for(int i = 0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_sync(thread_count, str, k, large_wind_kmer_am, start, minimized_h);
            }
        }
        
        // each thread gets 1 lwind
        thread_count = 86;
        for(int i =0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_sync(thread_count, str, k, large_wind_kmer_am, start, minimized_h);
            }
        }
        
        // some threads get 2 kmers, the rest get 1
        thread_count = 50;
        for(int i = 0; i < 10; i++){
            for(int i = 0; i < 4; i += 2){
                std::string str = test_strs[i].substr(start, 99);
                test_thread_sync(thread_count, str, k, large_wind_kmer_am, start, minimized_h);
            }
        }
        
    }
    
    SECTION("Full Testing"){
        unsigned thread_count = 4;
        unsigned k = 4;
        const int32_t large_wind_kmer_am = 8;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // the string changes
        // thread_count changes
        // start changes
        for(int i =0; i < 4; i++){
            for(int i = 0; i < 4; i += 2){
                for(int j = 4; j <= 64; j += 4){
                    thread_count = j;
                    for(int l = 0; l <= 96; l += 13){
                        start = l;
                        test_thread_sync(thread_count, test_strs[i], k, large_wind_kmer_am, start, minimized_h);
                    }
                }
            }
        }
    }
}
