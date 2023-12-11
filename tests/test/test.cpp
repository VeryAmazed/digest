//#include "nthash.hpp"
#include <catch2/catch_test_macros.hpp>
#include "digest/mod_minimizer.hpp"
#include "digest/window_minimizer.hpp"
#include "digest/syncmer.hpp"
#include <cstdint>
#include <fstream>

std::vector<std::string> test_strs;
// used to be first value was 1, but now k must be >= 4
unsigned ks[] = {4, 4, 7, 8, 9, 16, 25, 64};

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

void base_constructor(digest::Digester& dig, std::string& str, unsigned k, size_t pos, digest::MinimizedHashType minimized_h){
	INFO("String is: " << str);
	INFO("K is: " << k);
	INFO("Pos is: " << dig.get_pos());
	
	CHECK(strcmp(str.c_str(), dig.get_sequence()) == 0);
	CHECK(str.size() == dig.get_len());
	CHECK(dig.get_k() == k);
	CHECK(dig.get_minimized_h() == minimized_h);
	if(k <= str.size()){
		nthash::NtHash tHash(str, 1, k, pos);
		CHECK(dig.get_is_valid_hash() == tHash.roll());
		if(dig.get_is_valid_hash()){
			CHECK(dig.get_pos() == tHash.get_pos());
			INFO("ntHash pos is: " << tHash.get_pos());
			CHECK(dig.get_fhash() == tHash.get_forward_hash());
			CHECK(dig.get_rhash() == tHash.get_reverse_hash());
		}
	}else{
		CHECK(dig.get_is_valid_hash() == false);
	}
}

void base_dig_comp(digest::Digester& dig1, digest::Digester& dig2){
	CHECK(strcmp(dig1.get_sequence(), dig2.get_sequence()) == 0);
	CHECK(dig1.get_len() == dig2.get_len());
	CHECK(dig1.get_k() == dig2.get_k());
	CHECK(dig1.get_minimized_h() == dig2.get_minimized_h());
	CHECK(dig1.get_is_valid_hash() == dig2.get_is_valid_hash());
	if(dig1.get_is_valid_hash()){
		CHECK(dig1.get_chash() == dig2.get_chash());
		CHECK(dig1.get_fhash() == dig2.get_fhash());
		CHECK(dig1.get_rhash() == dig2.get_rhash());
	}
}

void base_dig_roll(digest::Digester& dig1, digest::Digester& dig2){
	while(dig1.get_is_valid_hash()){
		dig1.roll_one();
		dig2.roll_one();
		CHECK(dig1.get_fhash() == dig2.get_fhash());
		CHECK(dig1.get_rhash() == dig2.get_rhash());
		CHECK(dig1.get_pos() == dig2.get_pos());
	}
	CHECK(dig1.get_is_valid_hash() == dig2.get_is_valid_hash());
}

void ModMin_constructor(digest::ModMin& dig, std::string& str, unsigned k, size_t pos, digest::MinimizedHashType minimized_h, uint64_t mod, uint64_t congruence){
	base_constructor(dig, str, k, pos, minimized_h);
	CHECK(dig.get_mod() ==  mod);
	CHECK(dig.get_congruence() == congruence);
}

template <int32_t large_wind_kmer_am>
void WindowMin_constructor(digest::WindowMin<large_wind_kmer_am>& dig, std::string& str, unsigned k, size_t pos, digest::MinimizedHashType minimized_h){
	base_constructor(dig, str, k, pos, minimized_h);
	CHECK(dig.get_large_wind_kmer_am() == large_wind_kmer_am);
	// CHECK(dig.get_st_index() == 0);
	CHECK(dig.get_st_size() == 0);
	CHECK(dig.get_is_minimized() == false);
}

void ModMin_dig_comp(digest::ModMin& dig1, digest::ModMin& dig2){
	base_dig_comp(dig1, dig2);
	CHECK(dig1.get_mod() ==  dig2.get_mod());
	CHECK(dig1.get_congruence() == dig2.get_congruence());
	base_dig_roll(dig1, dig2);
}

template <int32_t k>
void WindowMin_roll_minimizers_comp(digest::WindowMin<k>& dig1, digest::WindowMin<k>& dig2){
	std::vector<size_t> vec1;
	std::vector<size_t> vec2;
	dig1.roll_minimizer(1000, vec1);
	dig2.roll_minimizer(1000, vec2);
	REQUIRE(vec1.size() == vec2.size());
	for(size_t i = 0; i < vec1.size(); i++){
		CHECK(vec1[i] == vec2[i]);
	}
}

template <int32_t k>
void Syncmer_roll_minimizers_comp(digest::Syncmer<k>& dig1, digest::Syncmer<k>& dig2){
	std::vector<std::pair<size_t, size_t>> vec1;
	std::vector<std::pair<size_t, size_t>> vec2;
	dig1.roll_minimizer(1000, vec1);
	dig2.roll_minimizer(1000, vec2);
	REQUIRE(vec1.size() == vec2.size());
	for(size_t i = 0; i < vec1.size(); i++){
		CHECK(vec1[i] == vec2[i]);
	}
}

template <int32_t k>
void WindowMin_dig_comp(digest::WindowMin<k>& dig1, digest::WindowMin<k>& dig2){
	base_dig_comp(dig1, dig2);
	CHECK(dig1.get_large_wind_kmer_am() == dig2.get_large_wind_kmer_am());
	// CHECK(dig1.get_st_index() == dig2.get_st_index());
	CHECK(dig1.get_st_size() == dig2.get_st_size());
	CHECK(dig1.get_is_minimized() == dig2.get_is_minimized());
	// need to use this because I need to check, or at least get some indication, of whether the two seg trees are the same
	WindowMin_roll_minimizers_comp(dig1, dig2);
}

template <int32_t k, int32_t l>
void Syncmer_dig_comp(digest::Syncmer<k>& dig1, digest::Syncmer<l>& dig2){
	base_dig_comp(dig1, dig2);
	CHECK(dig1.get_large_wind_kmer_am() == dig2.get_large_wind_kmer_am());
	// CHECK(dig1.get_st_index() == dig2.get_st_index());
	CHECK(dig1.get_st_size() == dig2.get_st_size());
	CHECK(dig1.get_is_minimized() == dig2.get_is_minimized());
	// need to use this because I need to check, or at least get some indication, of whether the two seg trees are the same
	Syncmer_roll_minimizers_comp(dig1, dig2);
}

void roll_one(digest::Digester& dig, std::string& str, unsigned k){
	INFO(str);
	INFO(k);
	nthash::NtHash tHash(str, 1, k, 0);
	uint64_t true_fhash;
	uint64_t true_rhash;
	uint64_t dig_fhash;
	uint64_t dig_rhash;
	bool worked = tHash.roll();
	while ((worked = tHash.roll())){
		dig.roll_one();
		CHECK(dig.get_is_valid_hash() == worked);
		if(worked){
			CHECK(dig.get_pos() == tHash.get_pos());
			true_fhash = tHash.get_forward_hash();
			true_rhash = tHash.get_reverse_hash();
			dig_fhash= dig.get_fhash();
			dig_rhash = dig.get_rhash();
			CHECK(dig_fhash == true_fhash);
			CHECK(dig_rhash == true_rhash);
		}
	}
	dig.roll_one();
	CHECK(dig.get_is_valid_hash() == worked);
}

void ModMin_roll_minimizer(digest::ModMin& dig, std::string& str, unsigned k, digest::MinimizedHashType minimized_h, uint32_t prime){
	nthash::NtHash tHash(str, 1, k, 0);
	std::vector<size_t> positions;
	std::vector<uint32_t> hashes;
	while(tHash.roll()){
		uint32_t temp;
		if(minimized_h == digest::MinimizedHashType::CANON){
			temp = *(tHash.hashes());
		}else if(minimized_h == digest::MinimizedHashType::FORWARD){
			temp = tHash.get_forward_hash();
		}else{
			temp = tHash.get_reverse_hash();
		}
		if(temp % prime == 0){
			positions.push_back(tHash.get_pos());
			hashes.push_back(temp%prime);
		}
		
	}
	std::vector<size_t> dig_positions;
	dig.roll_minimizer(400, dig_positions);
	REQUIRE(positions.size() == dig_positions.size());
	for(size_t i = 0; i < positions.size(); i++){
		CHECK(dig_positions[i] == positions[i]);
	}
}

template <int32_t large_wind_kmer_am>
void WindowMin_roll_minimizer(digest::WindowMin<large_wind_kmer_am>& dig, std::string& str, unsigned k, digest::MinimizedHashType minimized_h){
	nthash::NtHash tHash(str, 1, k, 0);
	std::vector<std::pair<uint32_t, size_t>> hashes;
	while(tHash.roll()){
		uint32_t temp;
		if(minimized_h == digest::MinimizedHashType::CANON){
			temp = *(tHash.hashes());
		}else if(minimized_h == digest::MinimizedHashType::FORWARD){
			temp = tHash.get_forward_hash();
		}else{
			temp = tHash.get_reverse_hash();
		}
		hashes.push_back(std::make_pair(temp, tHash.get_pos()));
	}

	std::vector<std::pair<uint32_t, size_t>> answers;
	std::pair<uint32_t, size_t> prev;
	for(size_t i =0; i + large_wind_kmer_am <= hashes.size(); i++){
		std::pair<uint32_t, size_t> temp_pair = hashes[i];
		for(uint j = 1; j < large_wind_kmer_am; j++){
			std::pair<uint32_t, size_t> curr = hashes[i+j];
			if(curr.first < temp_pair.first){
				temp_pair = curr;
			}else if(curr.first == temp_pair.first){
				if(curr.second > temp_pair.second){
					temp_pair = curr;
				}
			}
		}
		if(i == 0){
			prev = temp_pair;
			answers.push_back(temp_pair);
		}else{
			if(prev != temp_pair){
				prev = temp_pair;
				answers.push_back(temp_pair);
			}
		}
	}

	std::vector<size_t> wind_mins;
	dig.roll_minimizer(1000, wind_mins);
	REQUIRE(answers.size() == wind_mins.size());
	for(size_t i = 0; i < answers.size(); i++){
		CHECK(wind_mins[i] == answers[i].second);
	}
}

template <int32_t large_wind_kmer_am>
void Syncmer_roll_minimizer(digest::Syncmer<large_wind_kmer_am>& dig, std::string& str, unsigned k, digest::MinimizedHashType minimized_h){
	nthash::NtHash tHash(str, 1, k, 0);
	std::vector<std::pair<uint32_t, size_t>> hashes;
	while(tHash.roll()){
		uint32_t temp;
		if(minimized_h == digest::MinimizedHashType::CANON){
			temp = *(tHash.hashes());
		}else if(minimized_h == digest::MinimizedHashType::FORWARD){
			temp = tHash.get_forward_hash();
		}else{
			temp = tHash.get_reverse_hash();
		}
		hashes.push_back(std::make_pair(temp, tHash.get_pos()));
	}

	std::vector<std::pair<size_t, size_t>> answers;
	for(size_t i =0; i + large_wind_kmer_am <= hashes.size(); i++){
		uint32_t minAm = hashes[i].first;

		for(uint j = 1; j < large_wind_kmer_am; j++){
			minAm = std::min(minAm, hashes[i+j].first);
		}

		if(minAm == hashes[i].first || minAm == hashes[i+large_wind_kmer_am-1].first){
			answers.push_back(std::make_pair(hashes[i].second, hashes[i+large_wind_kmer_am-1].second));
		}
	}

	std::vector<std::pair<size_t, size_t>> wind_mins;
	dig.roll_minimizer(1000, wind_mins);
	REQUIRE(answers.size() == wind_mins.size());
	for(size_t i = 0; i < answers.size(); i++){
		CHECK(wind_mins[i].first == answers[i].first);
		CHECK(wind_mins[i].second == answers[i].second);
	}
}

void append_seq_compare(std::string& str1, std::string& str2, digest::Digester& dig, unsigned  k){
	INFO(str1);
	INFO(str2);
	INFO(str1.size());
	INFO(str2.size());
	INFO(k);
	
	std::string str3 = str1 + str2;
	nthash::NtHash tHash(str3, 1, k);
	std::vector<uint64_t> vec1;
	std::vector<size_t> positions1;
	while(tHash.roll()){
		vec1.push_back(*(tHash.hashes()));
		positions1.push_back(tHash.get_pos());
	}
	std::vector<uint64_t> vec2;
	std::vector<size_t> positions2;
	if(dig.get_is_valid_hash()){
		vec2.push_back(dig.get_chash());
		positions2.push_back(dig.get_pos());
		while(dig.roll_one()){
			vec2.push_back(dig.get_chash());
			positions2.push_back(dig.get_pos());
		}
	}
	dig.append_seq(str2);
	if(dig.get_is_valid_hash()){
		vec2.push_back(dig.get_chash());
		positions2.push_back(dig.get_pos());
		while(dig.roll_one()){
			vec2.push_back(dig.get_chash());
			positions2.push_back(dig.get_pos());
		}
	}
	REQUIRE(vec1.size() == vec2.size());
	for(size_t i =0; i < vec1.size(); i++){
		INFO(i);
		CHECK(vec1[i] == vec2[i]);
		CHECK(positions1[i] == positions2[i]);
	}
}

void append_seq_compare3(std::string& str1, std::string& str2, std::string str3, digest::Digester& dig, unsigned  k){
	INFO(str1);
	INFO(str2);
	INFO(str3);
	INFO(k);
	// Make sure to check positions too
	std::string str4 = str1 + str2 + str3;
	nthash::NtHash tHash(str4, 1, k);
	std::vector<uint64_t> vec1;
	std::vector<size_t> positions1;
	while(tHash.roll()){
		vec1.push_back(*(tHash.hashes()));
		positions1.push_back(tHash.get_pos());
	}
	std::vector<uint64_t> vec2;
	std::vector<size_t> positions2;
	if(dig.get_is_valid_hash()){
		vec2.push_back(dig.get_chash());
		positions2.push_back(dig.get_pos());
		while(dig.roll_one()){
			vec2.push_back(dig.get_chash());
			positions2.push_back(dig.get_pos());
		}
	}
	dig.append_seq(str2);
	if(dig.get_is_valid_hash()){
		vec2.push_back(dig.get_chash());
		positions2.push_back(dig.get_pos());
		while(dig.roll_one()){
			vec2.push_back(dig.get_chash());
			positions2.push_back(dig.get_pos());
		}
	}
	dig.append_seq(str3);
	if(dig.get_is_valid_hash()){
		vec2.push_back(dig.get_chash());
		positions2.push_back(dig.get_pos());
		while(dig.roll_one()){
			vec2.push_back(dig.get_chash());
			positions2.push_back(dig.get_pos());
		}
	}
	REQUIRE(vec1.size() == vec2.size());
	for(size_t i =0; i < vec1.size(); i++){
		INFO(i);
		CHECK(vec1[i] == vec2[i]);
		CHECK(positions1[i] == positions2[i]);
	}
}

void append_seq_small_cases(){
	std::string str1 = "CCGTGT";
	std::string str2 = "CCGNGT";
	std::string str3 = "AGCCTT";
	std::string str4 = "ANCCTT";
	std::string str5 = "A";

	digest::Digester* dig = new digest::ModMin(str1, 4, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare(str1, str3, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str2, 4, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare(str2, str4, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str2, 4, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare(str2, str3, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str2, 4, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare(str2, str5, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str1, 4, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare(str1, str5, *dig, 4);
	delete dig;
}

void append_seq_small_cases2(){
	std::string str1_good = "CATACCGGT";
	std::string str1_short = "TAG";
	std::string str1_badCh = "CATACNCGGT";

	std::string str2_good = "GTTCTCGCTT";
	std::string str2_badCh = "GTNTCTCGCTT";
	std::string str2A = "A";
	std::string str2_short = "TGGA";

	std::string str3_good = "CAACGACCGC";
	std::string str3_badCh = "NCAACGACCGC";

	digest::Digester* dig = new digest::ModMin(str1_good, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_good, str2_good, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_good, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_good, str2_badCh, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_good, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_good, str2A, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_short, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_short, str2A, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_badCh, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_badCh, str2A, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_good, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_good, str2_short, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_short, 6, 17, 0, 0, digest::MinimizedHashType::CANON);
	append_seq_compare3(str1_short, str2A, str3_badCh, *dig, 6);
	delete dig;
}
/*
	consider re-organizing this so this only tests the UM_Digester specific stuff
	like the constructor and roll_minimizer, but put the more general stuff, append_seq
	and roll_one in the general testing group
*/
TEST_CASE("Digester Testing"){
	setupStrings();
	// These use the ModMinimizer Class because Digester can't be instantiated, but correctness
	// doesn't depend on any of the fields of functions in the ModMinimizer Class
	SECTION("Base Constructor Special Cases"){
		unsigned k;
		digest::MinimizedHashType minimized_h;
		uint32_t mod, congruence;
		size_t pos;
		std::string str;
		// string is length 1, k = 1
		
		
		str = "AAAA";
		k = 4;
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = static_cast<digest::MinimizedHashType>(i);
			mod = 2;
			congruence = 1;
			digest::ModMin* dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h);
			ModMin_constructor(*dig, str, k, pos, minimized_h, mod, congruence);
			delete dig;
		}
		

		// string is length 1, k = 4
		str = "A";
		k = ks[1];
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = static_cast<digest::MinimizedHashType>(i);
			mod = 2;
			congruence = 1;
			digest::ModMin* dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h);
			ModMin_constructor(*dig, str, k, pos, minimized_h, mod, congruence);
			delete dig;
		}

		for(uint i =0; i < test_strs.size(); i++){
			for(int j =0; j < 8; j++){
				k = ks[j];
				for(int l =0; l < 16; l++){
					pos = l;
					for(int p =0; p < 3; p++){
						minimized_h = static_cast<digest::MinimizedHashType>(p);
						mod = 1e9+7;
						congruence = 0;
						digest::ModMin* dig = new digest::ModMin(test_strs[i], k, mod, congruence, pos, minimized_h);
						ModMin_constructor(*dig, test_strs[i], k, pos, minimized_h, mod, congruence);
						delete dig;
					}
				}
				
			}
		}

		// Throwing Exceptions
		// Shouldn't/Doesn't leak any memory
		// https://stackoverflow.com/questions/147572/will-the-below-code-cause-memory-leak-in-c
		
		str = "ACTGACTG";
		k = 4;
		pos = 0;
		minimized_h = digest::MinimizedHashType::CANON;
		mod = 1e9+7;
		congruence = 0;

		k = 0;
		digest::ModMin* dig;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		k = 4;
		
		// pos >= seq.size()
		pos = 8;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		pos = 0;

		/*
		// minimized_h > 2
		minimized_h = 3;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		minimized_h = 0;
		*/	
	}
	
	SECTION("Testing roll_one"){
		for(int i =0; i < 7; i++){
			for(int j =0; j < 8; j++){
				digest::ModMin* dig = new digest::ModMin(test_strs[i], ks[j], 1e9+7, 0, 0, digest::MinimizedHashType::FORWARD);
				roll_one(*dig, test_strs[i], ks[j]);
				delete dig;
			}
		}
	}

	SECTION("Testing append_seq()"){
		append_seq_small_cases();

		// Throws NotRolledTillEndException()
		digest::ModMin* dig = new digest::ModMin(test_strs[0], 4, 17);
		CHECK_THROWS_AS(dig->append_seq(test_strs[0]), digest::NotRolledTillEndException);
		delete dig;

		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					std::string str1 = test_strs[i].substr(0, l);
					std::string str2 = test_strs[i].substr(l, 100);
					digest::ModMin* dig = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, digest::MinimizedHashType::FORWARD);
					append_seq_compare(str1, str2, *dig, ks[j]);
					delete dig;
				}
			}
		}
		append_seq_small_cases2();
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					for(int r = 12; r < 85; r += 24){
						std::string str1 = test_strs[i].substr(0, l);
						std::string str2 = test_strs[i].substr(l, r);
						std::string str3 = test_strs[i].substr(l+r, 75);
						digest::ModMin* dig = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, digest::MinimizedHashType::FORWARD);
						append_seq_compare3(str1, str2, str3, *dig, ks[j]);
						delete dig;
					}
					
				}
			}
		}	
	}

	
	SECTION("Testing new_seq()"){
		unsigned k;
		std::string str;
		// string is length 1, k = 4
		str = "A";
		k = ks[1];
			
		digest::ModMin* dig1 = new digest::ModMin(test_strs[0], k, 1e9+7, 0, 0, digest::MinimizedHashType::CANON);
		dig1->new_seq(str, 0);
		base_constructor(*dig1, str, k, 0, digest::MinimizedHashType::CANON);
		delete dig1;
		
		// Throw BadConstructionException()
		dig1 = new digest::ModMin(test_strs[0], k, 1e9+7, 0, 0, digest::MinimizedHashType::CANON);
		CHECK_THROWS_AS(dig1->new_seq(test_strs[0], 500), digest::BadConstructionException);
		delete dig1;

		for(uint i =0; i < test_strs.size(); i += 2){
			for(int j = 0; j < 32; j += 8){
				digest::ModMin* dig = new digest::ModMin(test_strs[1], ks[3], 1e9+7, 0, 0, digest::MinimizedHashType::CANON);
				dig->new_seq(test_strs[i], j);
				base_constructor(*dig, test_strs[i], ks[3], j, digest::MinimizedHashType::CANON);
				delete dig;
			}
		}

		for(uint i =0; i < test_strs.size(); i += 2){
			for(int l = 13; l <= 78; l +=13 ){
				digest::ModMin* dig = new digest::ModMin(test_strs[5], ks[3], 1e9+7, 0, 0, digest::MinimizedHashType::CANON);
				int ind = 0;
				while(ind < l && dig->roll_one()){
					ind++;
				}
				dig->new_seq(test_strs[i], 0);
				base_constructor(*dig, test_strs[i], ks[3], 0, digest::MinimizedHashType::CANON);
				delete dig;
			}
			
		}
		
		// new_seq when deque has stuff in it
		dig1 = new digest::ModMin(test_strs[2], 8, 17, 0, 0, digest::MinimizedHashType::CANON);
		std::vector<size_t> vec;
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->append_seq(test_strs[2]);
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->new_seq(test_strs[4], 0);
		base_constructor(*dig1, test_strs[4], 8, 0, digest::MinimizedHashType::CANON);
		delete dig1;

		// new_seq when deque has stuff in it and a new hash can't be properly initialized
		std::string bad_str = "TTACTNGTACCTG";
		dig1 = new digest::ModMin(test_strs[2], 8, 17, 0, 0, digest::MinimizedHashType::CANON);
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->append_seq(test_strs[2]);
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->new_seq(bad_str, 0);
		base_constructor(*dig1, bad_str, 8, 0, digest::MinimizedHashType::CANON);
		delete dig1;
	}
}

TEST_CASE("ModMin Testing"){
	setupStrings();

	SECTION("Testing Constructors"){
		unsigned k;
		digest::MinimizedHashType minimized_h;
		uint32_t mod, congruence;
		size_t pos;
		std::string str;

		// Throwing Exceptions
		// Shouldn't/Doesn't leak any memory
		// https://stackoverflow.com/questions/147572/will-the-below-code-cause-memory-leak-in-c
		
		
		str = "ACTGACTG";
		k = 4;
		pos = 0;
		minimized_h = digest::MinimizedHashType::CANON;
		digest::ModMin* dig;

		// mod >= congruence
		mod = 2;
		congruence = 2;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadModException);
		
	}
	
	// maybe move this into an entirely new test case, and make this big thing just tests for the Dig class
	SECTION("Testing roll_minimizer(). The one that takes no parameters"){
		uint32_t prime = 17;
		for(int i =0; i < 7; i += 2){
			for(int j =0; j < 8; j++){
				for(int l = 0; l < 3; l++){
					digest::ModMin* dig = new digest::ModMin(test_strs[i], ks[j], prime, 0, 0, static_cast<digest::MinimizedHashType>(l));
					ModMin_roll_minimizer(*dig, test_strs[i], ks[j], static_cast<digest::MinimizedHashType>(l), prime);
					delete dig;
				}
				
			}
		}
	}
	
	SECTION("Testing Copy Constructor"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					digest::ModMin* dig1 = new digest::ModMin(test_strs[i], ks[j], 1e9+7, 0, l, digest::MinimizedHashType::FORWARD);
					digest::ModMin* dig2 = new digest::ModMin(*dig1);
					ModMin_dig_comp(*dig1, *dig2);
					delete dig1;
					delete dig2;
				}
			}
		}

		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					std::string str1 = test_strs[i].substr(0, l);
					std::string str2 = test_strs[i].substr(l, 100);
					digest::ModMin* dig1 = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, digest::MinimizedHashType::FORWARD);
					std::vector<size_t> vec;
					dig1->roll_minimizer(1000, vec);
					dig1->append_seq(str2);
					digest::ModMin* dig2 = new digest::ModMin(*dig1);
					ModMin_dig_comp(*dig1, *dig2);
					delete dig1;
					delete dig2;
				}
			}
		}
	}

	SECTION("Testing Assignment Operator"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					digest::ModMin* dig1 = new digest::ModMin(test_strs[i], ks[j], 1e9+7, 0, l, digest::MinimizedHashType::FORWARD);
					digest::ModMin* dig2 = new digest::ModMin(test_strs[1], 99, 98765, 3, 0, digest::MinimizedHashType::REVERSE);
					*dig2 = *dig1;
					ModMin_dig_comp(*dig1, *dig2);
					delete dig1;
					delete dig2;
				}
			}
		}

		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					std::string str1 = test_strs[i].substr(0, l);
					std::string str2 = test_strs[i].substr(l, 100);
					digest::ModMin* dig1 = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, digest::MinimizedHashType::FORWARD);
					std::vector<size_t> vec;
					dig1->roll_minimizer(1000, vec);
					dig1->append_seq(str2);
					digest::ModMin* dig2 = new digest::ModMin(test_strs[1], 99, 98765, 3, 0, digest::MinimizedHashType::REVERSE);
					*dig2 = *dig1;
					ModMin_dig_comp(*dig1, *dig2);
					delete dig1;
					delete dig2;
				}
			}
		}
	}
}

TEST_CASE("WindowMin Testing"){
	SECTION("Constructor Testing"){
		unsigned k;
		digest::MinimizedHashType minimized_h;
		size_t pos;
		std::string str;

		// Throwing Exceptions
		// Shouldn't/Doesn't leak any memory
		// https://stackoverflow.com/questions/147572/will-the-below-code-cause-memory-leak-in-c
		
		str = "ACTGACTG";
		k = 4;
		pos = 0;
		minimized_h = digest::MinimizedHashType::CANON;
		digest::WindowMin<0>* dig1;
		CHECK_THROWS_AS((dig1 = new digest::WindowMin<0>(str, k, pos, minimized_h)), digest::BadWindowSizeException);

		for(uint i =0; i < test_strs.size(); i++){
			k = 4;
			pos = 0;
			minimized_h = digest::MinimizedHashType::CANON;

			# define TEST_CONSTRUCTOR_0(j) \
				digest::WindowMin<j>* dig = new digest::WindowMin<j>(test_strs[i], k, pos, minimized_h); \
				WindowMin_constructor(*dig, test_strs[i], k, pos, minimized_h); \
				delete dig;

			{ TEST_CONSTRUCTOR_0(1) }
			{ TEST_CONSTRUCTOR_0(2) }
			{ TEST_CONSTRUCTOR_0(3) }
			{ TEST_CONSTRUCTOR_0(4) }
			{ TEST_CONSTRUCTOR_0(13) }
			{ TEST_CONSTRUCTOR_0(14) }
			{ TEST_CONSTRUCTOR_0(18) }
			{ TEST_CONSTRUCTOR_0(20) }
			{ TEST_CONSTRUCTOR_0(31) }
			{ TEST_CONSTRUCTOR_0(32) }
			{ TEST_CONSTRUCTOR_0(63) }
			{ TEST_CONSTRUCTOR_0(64) }
		}
	}

	SECTION("roll_minimizer() testing"){
		for(int i =0; i < 7; i += 2){
			//std::cout << test_strs[i] << std::endl;
			for(int j =0; j < 8; j++){
				for(int l = 0; l < 3; l++){
					# define TEST_ROLL_0(m) \
						digest::WindowMin<m>* dig = new digest::WindowMin<m>(test_strs[i], ks[j], 0, static_cast<digest::MinimizedHashType>(l)); \
						WindowMin_roll_minimizer(*dig, test_strs[i], ks[j], static_cast<digest::MinimizedHashType>(l)); \
						delete dig;
					{ TEST_ROLL_0(1) }
					{ TEST_ROLL_0(2) }
					{ TEST_ROLL_0(3) }
					{ TEST_ROLL_0(4) }
					{ TEST_ROLL_0(13) }
					{ TEST_ROLL_0(14) }
					{ TEST_ROLL_0(18) }
					{ TEST_ROLL_0(20) }
					{ TEST_ROLL_0(31) }
					{ TEST_ROLL_0(32) }
				}
			}
		}
	}
	/*
		the below also inadverntently tests how append_seq (only the case that there are 2 sequences involved total) 
		works with roll_minimizer for WindowMin. In theory this shouldn't be needed and also can't be considered "thorough", but it is extra assurance.
	*/
	SECTION("Testing Copy Constructor"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_COPY_0(k) \
						digest::WindowMin<k> dig1(test_strs[i], ks[j], l, digest::MinimizedHashType::FORWARD); \
						digest::WindowMin<k> dig2(dig1); \
						WindowMin_dig_comp(dig1, dig2);
					{ TEST_COPY_0(1) }
					{ TEST_COPY_0(2) }
					{ TEST_COPY_0(3) }
					{ TEST_COPY_0(4) }
					{ TEST_COPY_0(13) }
					{ TEST_COPY_0(14) }
					{ TEST_COPY_0(18) }
					{ TEST_COPY_0(20) }
					{ TEST_COPY_0(31) }
					{ TEST_COPY_0(32) }
				}
			}
		}

		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_COPY_1(k) \
						std::string str1 = test_strs[i].substr(0, l); \
						std::string str2 = test_strs[i].substr(l, 100); \
						digest::WindowMin<k> dig1(str1, ks[j], 0, digest::MinimizedHashType::FORWARD); \
						std::vector<size_t> vec; \
						dig1.roll_minimizer(1000, vec); \
						dig1.append_seq(str2); \
						digest::WindowMin<k> dig2(dig1); \
						WindowMin_dig_comp(dig1, dig2);
					{ TEST_COPY_1(1) }
					{ TEST_COPY_1(2) }
					{ TEST_COPY_1(3) }
					{ TEST_COPY_1(4) }
					{ TEST_COPY_1(13) }
					{ TEST_COPY_1(14) }
					{ TEST_COPY_1(18) }
					{ TEST_COPY_1(20) }
					{ TEST_COPY_1(31) }
					{ TEST_COPY_1(32) }
				}
			}
		}
	}

	SECTION("Testing Assignment Operator"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_COPY_2(k) \
						digest::WindowMin<k> dig1(test_strs[i], ks[j], l, digest::MinimizedHashType::FORWARD); \
						digest::WindowMin<k> dig2(test_strs[1], 99, 0, digest::MinimizedHashType::REVERSE); \
						dig2 = dig1; \
						WindowMin_dig_comp(dig1, dig2);
					{ TEST_COPY_2(1) }
					{ TEST_COPY_2(2) }
					{ TEST_COPY_2(3) }
					{ TEST_COPY_2(4) }
					{ TEST_COPY_2(13) }
					{ TEST_COPY_2(14) }
					{ TEST_COPY_2(18) }
					{ TEST_COPY_2(20) }
					{ TEST_COPY_2(31) }
					{ TEST_COPY_2(32) }
				}
			}
		}


		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_COPY_3(m) \
						std::string str1 = test_strs[i].substr(0, l); \
						std::string str2 = test_strs[i].substr(l, 100); \
						digest::WindowMin<m> dig1(str1, ks[j], 0, digest::MinimizedHashType::FORWARD); \
						std::vector<size_t> vec; \
						dig1.roll_minimizer(1000, vec); \
						dig1.append_seq(str2); \
						digest::WindowMin<m> dig2(test_strs[1], 35, 0, digest::MinimizedHashType::REVERSE); \
						dig2 = dig1; \
						WindowMin_dig_comp(dig1, dig2);
					{ TEST_COPY_3(1) }
					{ TEST_COPY_3(2) }
					{ TEST_COPY_3(3) }
					{ TEST_COPY_3(4) }
					{ TEST_COPY_3(13) }
					{ TEST_COPY_3(14) }
					{ TEST_COPY_3(18) }
					{ TEST_COPY_3(20) }
					{ TEST_COPY_3(31) }
					{ TEST_COPY_3(32) }
				}
			}
		}
	}
}

TEST_CASE("Syncmer Testing"){
	// Syncmer and WindowMinimizers have all the same class members so I can just use the WindowMin tests for Constructor and be ok
	SECTION("Constructor Testing"){
		
		unsigned k;
		digest::MinimizedHashType minimized_h;
		size_t pos;
		std::string str;

		// Throwing Exceptions
		// Shouldn't/Doesn't leak any memory
		// https://stackoverflow.com/questions/147572/will-the-below-code-cause-memory-leak-in-c
		
		str = "ACTGACTG";
		k = 4;
		pos = 0;
		minimized_h = digest::MinimizedHashType::CANON;
		digest::Syncmer<0>* dig1;
		CHECK_THROWS_AS((dig1 = new digest::Syncmer<0>(str, k, pos, minimized_h)), digest::BadWindowSizeException);

		for(uint i =0; i < test_strs.size(); i++){
			# define TEST_SYNCON(j) \
				k = 4; \
				pos = 0; \
				minimized_h = digest::MinimizedHashType::CANON; \
				 \
				digest::Syncmer<j> dig(test_strs[i], k, pos, minimized_h); \
				WindowMin_constructor(dig, test_strs[i], k, pos, minimized_h);

			{ TEST_SYNCON(1) }
			{ TEST_SYNCON(2) }
			{ TEST_SYNCON(3) }
			{ TEST_SYNCON(4) }
			{ TEST_SYNCON(13) }
			{ TEST_SYNCON(14) }
			{ TEST_SYNCON(18) }
			{ TEST_SYNCON(20) }
			{ TEST_SYNCON(31) }
			{ TEST_SYNCON(32) }
			{ TEST_SYNCON(63) }
			{ TEST_SYNCON(64) }
		}
	}

	SECTION("roll_minimizer() testing"){
		for(int i =0; i < 7; i += 2){
			//std::cout << test_strs[i] << std::endl;
			for(int j =0; j < 8; j++){
				for(int l = 0; l < 3; l++){
					;
					# define TEST_SYNCROLL(m) \
						digest::Syncmer<m> dig(test_strs[i], ks[j], 0, static_cast<digest::MinimizedHashType>(l)); \
						Syncmer_roll_minimizer(dig, test_strs[i], ks[j], static_cast<digest::MinimizedHashType>(l)); \

					{ TEST_SYNCROLL(1) }
					{ TEST_SYNCROLL(2) }
					{ TEST_SYNCROLL(3) }
					{ TEST_SYNCROLL(4) }
					{ TEST_SYNCROLL(13) }
					{ TEST_SYNCROLL(14) }
					{ TEST_SYNCROLL(18) }
					{ TEST_SYNCROLL(20) }
					{ TEST_SYNCROLL(31) }
					{ TEST_SYNCROLL(32) }
				}
			}

		}
	}

	/*
		the below also inadverntently tests how append_seq (only the case that there are 2 sequences involved total) 
		works with roll_minimizer for WindowMin. In theory this shouldn't be needed and also can't be considered "thorough", but it is extra assurance.
	*/
	SECTION("Testing Copy Constructor"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_SYNCOPY_0(m) \
						digest::Syncmer<m> dig1(test_strs[i], ks[j], l, digest::MinimizedHashType::FORWARD); \
						digest::Syncmer<m> dig2(dig1); \
						Syncmer_dig_comp(dig1, dig2);

					{ TEST_SYNCOPY_0(1) }
					{ TEST_SYNCOPY_0(2) }
					{ TEST_SYNCOPY_0(3) }
					{ TEST_SYNCOPY_0(4) }
					{ TEST_SYNCOPY_0(13) }
					{ TEST_SYNCOPY_0(14) }
					{ TEST_SYNCOPY_0(18) }
					{ TEST_SYNCOPY_0(20) }
					{ TEST_SYNCOPY_0(31) }
					{ TEST_SYNCOPY_0(32) }
				}
			}
		}

		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_SYNCOPY_1(m) \
						std::string str1 = test_strs[i].substr(0, l); \
						std::string str2 = test_strs[i].substr(l, 100); \
						digest::Syncmer<m> dig1(str1, ks[j], 0, digest::MinimizedHashType::FORWARD); \
						std::vector<size_t> vec; \
						dig1.roll_minimizer(1000, vec); \
						dig1.append_seq(str2); \
						digest::Syncmer<m> dig2(dig1); \
						Syncmer_dig_comp(dig1, dig2);

					{ TEST_SYNCOPY_1(1) }
					{ TEST_SYNCOPY_1(2) }
					{ TEST_SYNCOPY_1(3) }
					{ TEST_SYNCOPY_1(4) }
					{ TEST_SYNCOPY_1(13) }
					{ TEST_SYNCOPY_1(14) }
					{ TEST_SYNCOPY_1(18) }
					{ TEST_SYNCOPY_1(20) }
					{ TEST_SYNCOPY_1(31) }
					{ TEST_SYNCOPY_1(32) }
					
				}
			}
		}
	}

	// had to modify !!!
	SECTION("Testing Assignment Operator"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_ASSIGNMENT_0(m) \
						digest::Syncmer<m>* dig1 = new digest::Syncmer<m>(test_strs[i], ks[j],l, digest::MinimizedHashType::FORWARD);\
						digest::Syncmer<m>* dig2 = new digest::Syncmer<m>(test_strs[1], 99, 0, digest::MinimizedHashType::REVERSE);\
						*dig2 = *dig1;\
						Syncmer_dig_comp(*dig1, *dig2);\
						delete dig1;\
						delete dig2;

					{ TEST_ASSIGNMENT_0(1) }
					{ TEST_ASSIGNMENT_0(2) }
					{ TEST_ASSIGNMENT_0(3) }
					{ TEST_ASSIGNMENT_0(4) }
					{ TEST_ASSIGNMENT_0(13) }
					{ TEST_ASSIGNMENT_0(14) }
					{ TEST_ASSIGNMENT_0(18) }
					{ TEST_ASSIGNMENT_0(20) }
					{ TEST_ASSIGNMENT_0(31) }
					{ TEST_ASSIGNMENT_0(32) }
					
				}
			}
		}

		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					# define TEST_ASSIGNMENT_1(m) \
						std::string str1 = test_strs[i].substr(0, l); \
						std::string str2 = test_strs[i].substr(l, 100); \
						digest::Syncmer<m>* dig1 = new digest::Syncmer<m>(str1, ks[j], 0, digest::MinimizedHashType::FORWARD); \
						std::vector<size_t> vec; \
						dig1->roll_minimizer(1000, vec); \
						dig1->append_seq(str2); \
						digest::Syncmer<m>* dig2 = new digest::Syncmer<m>(test_strs[1], 35, 0, digest::MinimizedHashType::REVERSE); \
						*dig2 = *dig1; \
						Syncmer_dig_comp(*dig1, *dig2); \
						delete dig1; \
						delete dig2;
					{ TEST_ASSIGNMENT_1(1) }
					{ TEST_ASSIGNMENT_1(2) }
					{ TEST_ASSIGNMENT_1(3) }
					{ TEST_ASSIGNMENT_1(4) }
					{ TEST_ASSIGNMENT_1(13) }
					{ TEST_ASSIGNMENT_1(14) }
					{ TEST_ASSIGNMENT_1(18) }
					{ TEST_ASSIGNMENT_1(20) }
					{ TEST_ASSIGNMENT_1(31) }
					{ TEST_ASSIGNMENT_1(32) }
				}
			}
		}
	}
	
}
//
// #include <iostream>
//
// template <int k>
// class MyClass {
// public:
//     void display() {
//         std::cout << "Value of k: " << k << std::endl;
//     }
// };
//
// // Template specialization to handle the base case
// template <>
// class MyClass<0> {
// public:
//     void display() {
//         // Do nothing or handle the base case as needed
//     }
// };
//
// // Recursive template to instantiate objects with values from 1 to 32
// template <int i>
// struct InstantiateObjects {
//     static void instantiate() {
//         MyClass<i> obj;
//         obj.display();
//         InstantiateObjects<i - 1>::instantiate();
//     }
// };
//
// // Template specialization to handle the base case of recursion
// template <>
// struct InstantiateObjects<0> {
//     static void instantiate() {
//         // Do nothing or handle the base case as needed
//     }
// };
//
// int main() {
//     // Instantiate objects with values from 1 to 32
//     InstantiateObjects<32>::instantiate();
//
//     return 0;
// }
//
