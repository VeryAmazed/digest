#include "nthash.hpp"
#include <catch2/catch_test_macros.hpp>
#include "mod_minimizer.hpp"
#include <fstream>

std::vector<std::string> test_strs;
unsigned ks[] = {1, 4, 7, 8, 9, 16, 25, 64};

void setupStrings(){
	std::string str;
	std::fstream fs;

	fs.open("../test_strings/A.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/a_lowercase.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/salmonella_enterica.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/salmonella_lowercase.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/random.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/random_lowercase.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/N.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();
}

void base_constructor(digest::Digester& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h){
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

void ModMin_constructor(digest::ModMin& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h, uint64_t mod, uint64_t congruence){
	base_constructor(dig, str, k, pos, minimized_h);
	CHECK(dig.get_mod() ==  mod);
	CHECK(dig.get_congruence() == congruence);
}

void ModMin_dig_comp(digest::ModMin& dig1, digest::ModMin& dig2){
	base_dig_comp(dig1, dig2);
	CHECK(dig1.get_mod() ==  dig2.get_mod());
	CHECK(dig1.get_congruence() == dig2.get_congruence());
	base_dig_roll(dig1, dig2);
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
	while(worked = tHash.roll()){
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

void ModMin_roll_minimizer(digest::Digester& dig, std::string& str, unsigned k, unsigned minimized_h, uint64_t prime){
	nthash::NtHash tHash(str, 1, k, 0);
	std::vector<size_t> positions;
	std::vector<uint64_t> hashes;
	while(tHash.roll()){
		uint64_t temp;
		if(minimized_h == 0){
			temp = *(tHash.hashes());
		}else if(minimized_h == 1){
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

void append_seq_compare(std::string& str1, std::string& str2, digest::Digester& dig, unsigned  k){
	INFO(str1);
	INFO(str2);
	INFO(str1.size());
	INFO(str2.size());
	INFO(k);
	// Make sure to check positions too
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

	digest::Digester* dig = new digest::ModMin(str1, 4, 17, 0, 0, 0);
	append_seq_compare(str1, str3, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str2, 4, 17, 0, 0, 0);
	append_seq_compare(str2, str4, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str2, 4, 17, 0, 0, 0);
	append_seq_compare(str2, str3, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str2, 4, 17, 0, 0, 0);
	append_seq_compare(str2, str5, *dig, 4);
	delete dig;

	dig = new digest::ModMin(str1, 4, 17, 0, 0, 0);
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

	digest::Digester* dig = new digest::ModMin(str1_good, 6, 17, 0, 0, 0);
	append_seq_compare3(str1_good, str2_good, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_good, 6, 17, 0, 0, 0);
	append_seq_compare3(str1_good, str2_badCh, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_good, 6, 17, 0, 0, 0);
	append_seq_compare3(str1_good, str2A, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_short, 6, 17, 0, 0, 0);
	append_seq_compare3(str1_short, str2A, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_badCh, 6, 17, 0, 0, 0);
	append_seq_compare3(str1_badCh, str2A, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_good, 6, 17, 0, 0, 0);
	append_seq_compare3(str1_good, str2_short, str3_good, *dig, 6);
	delete dig;

	dig = new digest::ModMin(str1_short, 6, 17, 0, 0, 0);
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
		unsigned k, minimized_h;
		uint64_t mod, congruence;
		size_t pos;
		std::string str;
		// string is length 1, k = 1
		str = "A";
		k = ks[0];
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = i;
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
			minimized_h = i;
			mod = 2;
			congruence = 1;
			digest::ModMin* dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h);
			ModMin_constructor(*dig, str, k, pos, minimized_h, mod, congruence);
			delete dig;
		}

		// Throwing Exceptions
		// Shouldn't/Doesn't leak any memory
		// https://stackoverflow.com/questions/147572/will-the-below-code-cause-memory-leak-in-c
		
		str = "ACTGACTG";
		k = 2;
		pos = 0;
		minimized_h = 0;
		mod = 1e9+7;
		congruence = 0;
		// k = 0
		k = 0;
		digest::ModMin* dig;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		k = 2;
		
		// pos >= seq.size()
		pos = 8;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		pos = 0;

		// minimized_h > 2
		minimized_h = 3;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		minimized_h = 0;	
	}
	
	SECTION("Testing roll_one"){
		for(int i =0; i < 7; i++){
			for(int j =0; j < 8; j++){
				digest::ModMin* dig = new digest::ModMin(test_strs[i], ks[j], 1e9+7, 0, 0, 1);
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
					digest::ModMin* dig = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, 1);
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
						digest::ModMin* dig = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, 1);
						append_seq_compare3(str1, str2, str3, *dig, ks[j]);
						delete dig;
					}
					
				}
			}
		}	
	}

	
	SECTION("Testing new_seq()"){
		unsigned k, minimized_h;
		std::string str;
		// string is length 1, k = 4
		str = "A";
		k = ks[1];
			
		digest::ModMin* dig1 = new digest::ModMin(test_strs[0], k, 1e9+7, 0, 0, 0);
		dig1->new_seq(str, 0);
		base_constructor(*dig1, str, k, 0, 0);
		delete dig1;
		
		// Throw BadConstructionException()
		dig1 = new digest::ModMin(test_strs[0], k, 1e9+7, 0, 0, 0);
		CHECK_THROWS_AS(dig1->new_seq(test_strs[0], 500), digest::BadConstructionException);
		delete dig1;

		for(int i =0; i < test_strs.size(); i += 2){
			for(int j = 0; j < 32; j += 8){
				digest::ModMin* dig = new digest::ModMin(test_strs[1], ks[3], 1e9+7, 0, 0, 0);
				dig->new_seq(test_strs[i], j);
				base_constructor(*dig, test_strs[i], ks[3], j, 0);
				delete dig;
			}
		}

		for(int i =0; i < test_strs.size(); i += 2){
			for(int l = 13; l <= 78; l +=13 ){
				digest::ModMin* dig = new digest::ModMin(test_strs[5], ks[3], 1e9+7, 0, 0, 0);
				int ind = 0;
				while(ind < l && dig->roll_one()){
					ind++;
				}
				dig->new_seq(test_strs[i], 0);
				base_constructor(*dig, test_strs[i], ks[3], 0, 0);
				delete dig;
			}
			
		}
		
		// new_seq when deque has stuff in it
		dig1 = new digest::ModMin(test_strs[2], 8, 17, 0, 0, 0);
		std::vector<size_t> vec;
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->append_seq(test_strs[2]);
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->new_seq(test_strs[4], 0);
		base_constructor(*dig1, test_strs[4], 8, 0, 0);
		delete dig1;

		// new_seq when deque has stuff in it and a new hash can't be properly initialized
		std::string bad_str = "TTACTNGTACCTG";
		dig1 = new digest::ModMin(test_strs[2], 8, 17, 0, 0, 0);
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->append_seq(test_strs[2]);
		dig1->roll_minimizer(1000, vec);
		vec.clear();
		dig1->new_seq(bad_str, 0);
		base_constructor(*dig1, bad_str, 8, 0, 0);
		delete dig1;
	}
}

TEST_CASE("ModMin Testing"){
	setupStrings();

	SECTION("Testing Constructors"){
		unsigned k, minimized_h;
		uint64_t mod, congruence;
		size_t pos;
		std::string str;
		// string is length 1, k = 1
		
		for(int i =0; i < test_strs.size(); i++){
			for(int j =0; j < 8; j++){
				k = ks[j];
				for(int l =0; l < 16; l++){
					pos = l;
					for(int p =0; p < 3; p++){
						minimized_h = p;
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
		k = 2;
		pos = 0;
		minimized_h = 0;
		digest::ModMin* dig;

		// mod >= congruence
		mod = 2;
		congruence = 2;
		CHECK_THROWS_AS(dig = new digest::ModMin(str, k, mod, congruence, pos, minimized_h), digest::BadModException);
		
	}
	
	// maybe move this into an entirely new test case, and make this big thing just tests for the Dig class
	SECTION("Testing roll_minimizer(). The one that takes no parameters"){
		uint64_t prime = 17;
		for(int i =0; i < 7; i += 2){
			for(int j =0; j < 8; j++){
				for(int l = 0; l < 3; l++){
					digest::ModMin* dig = new digest::ModMin(test_strs[i], ks[j], prime, 0, 0, l);
					ModMin_roll_minimizer(*dig, test_strs[i], ks[j], l, prime);
					delete dig;
				}
				
			}
		}
	}
	
	SECTION("Testing Copy Constructor"){
		for(int i =0; i < 7; i +=2){
			for(int j =0; j < 8; j++){
				for(int l = 15; l < 91; l += 15){
					digest::ModMin* dig1 = new digest::ModMin(test_strs[i], ks[j], 1e9+7, 0, l, 1);
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
					digest::ModMin* dig1 = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, 1);
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
					digest::ModMin* dig1 = new digest::ModMin(test_strs[i], ks[j], 1e9+7, 0, l, 1);
					digest::ModMin* dig2 = new digest::ModMin(test_strs[1], 99, 98765, 3, 0, 2);
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
					digest::ModMin* dig1 = new digest::ModMin(str1, ks[j], 1e9+7, 0, 0, 1);
					std::vector<size_t> vec;
					dig1->roll_minimizer(1000, vec);
					dig1->append_seq(str2);
					digest::ModMin* dig2 = new digest::ModMin(test_strs[1], 99, 98765, 3, 0, 2);
					*dig2 = *dig1;
					ModMin_dig_comp(*dig1, *dig2);
					delete dig1;
					delete dig2;
				}
			}
		}
	}
}
