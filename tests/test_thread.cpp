#include <catch2/catch_test_macros.hpp>
#include "digest/thread_out.hpp"
#include <fstream>
#include <iostream>

std::vector<std::string> test_strs;

void setupStrings(){
	std::string files[] = {
		"../tests/test_strings/A.txt",
		"../tests/test_strings/a_lowercase.txt",
		"../tests/test_strings/salmonella_enterica.txt",
		"../tests/test_strings/salmonella_lowercase.txt",
		"../tests/test_strings/random.txt",
		"../tests/test_strings/random_lowercase.txt",
		"../tests/test_strings/N.txt",
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

void test_thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    std::string str, unsigned k, uint64_t mod, uint64_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
    std::vector<size_t> single_thread;
    digest::ModMin dig(str, k, mod, congruence, start, minimized_h);
    dig.roll_minimizer(str.size(), single_thread);
    thread_out::thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h);
    std::vector<size_t> multi_thread = multi_to_single_vec(vec);
    std::cout << single_thread.size() << " " << multi_thread.size() << std::endl;
    for(size_t i = 0; i < single_thread.size(); i++){
		std::cout << single_thread[i] << " ";
	}
    std::cout << std::endl;
    for(size_t i = 0; i < multi_thread.size(); i++){
		std::cout << multi_thread[i] << " ";
	}
    std::cout << std::endl;
    REQUIRE(single_thread.size() == multi_thread.size());
	for(size_t i = 0; i < single_thread.size(); i++){
		CHECK(single_thread[i] == multi_thread[i]);
	}
}

TEST_CASE("thread_out function testing"){
    setupStrings();
    SECTION("Throw Errors"){
        std::string str = "ACTGACTG";
        unsigned thread_count = 4;
        std::vector<std::vector<size_t>> vec;
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
        std::vector<std::vector<size_t>> vec;
        unsigned k = 4;
        uint64_t mod = 3;
        uint64_t congruence = 0;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // only 1 thread
        thread_count = 1;
        for(int i = 0; i < 5; i += 2){
            std::string str = test_strs[i].substr(start, 99);
            test_thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h);
        }
        std::cout << "Pass 1 thread" << std::endl;
        // each thread gets 1 kmer
        thread_count = 96;
        for(int i = 0; i < 5; i += 2){
            std::string str = test_strs[i].substr(start, 99);
            test_thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h);
        }
        std::cout << "Pass 1 kmer per thread" << std::endl;
        // some threads get 2 kmers, the rest get 1
        thread_count = 50;
        for(int i = 0; i < 5; i += 2){
            std::string str = test_strs[i].substr(start, 99);
            test_thread_mod(thread_count, vec, str, k, mod, congruence, start, minimized_h);
        }
        std::cout << "pass some kmers get 2" << std::endl;
    }

    SECTION("Full Testing"){
        unsigned thread_count = 4;
        std::vector<std::vector<size_t>> vec;
        unsigned k = 4;
        uint64_t mod = 17;
        uint64_t congruence = 0;
        size_t start = 0;
        digest::MinimizedHashType minimized_h = digest::MinimizedHashType::CANON;
        // the string changes
        // thread_count changes
        // start changes
        for(int i = 0; i < 5; i += 2){
            for(int j = 4; j <= 64; j += 4){
                thread_count = j;
                for(int l = 0; l <= 96; l += 13){
                    start = l;
                    test_thread_mod(thread_count, vec, test_strs[i], k, mod, congruence, start, minimized_h);
                }
            }
        }
    }
}