#include "nthash.hpp"
#include <catch2/catch_test_macros.hpp>
#include "uni_mini.hpp"
#include <fstream>

std::vector<std::string> test_strs;
unsigned ks[] = {1, 4, 9, 89, 128};

void setupStrings(){
	std::string str;
	std::fstream fs;

	fs.open("../test_strings/A.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/N.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/a_lowercase.txt", std::fstream::in);
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

	fs.open("../test_strings/salmonella_enterica.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();

	fs.open("../test_strings/salmonella_lowercase.txt", std::fstream::in);
	fs >> str;
	test_strs.push_back(str);
	fs.close();
}

void base_constructor_stdstr(digest::Digester& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h){
	CHECK(strcmp(str.c_str(), dig.get_sequence()) == 0);
	CHECK(str.size() == dig.get_len());
	CHECK(dig.get_k() == k);
	CHECK(dig.get_pos() == pos);
	CHECK(dig.get_minimized_h() == minimized_h);
	CHECK(dig.get_rolled() == false);
	CHECK_THROWS_AS(dig.get_chash(), digest::NotRolledException);
	CHECK_THROWS_AS(dig.get_fhash(), digest::NotRolledException);
	CHECK_THROWS_AS(dig.get_rhash(), digest::NotRolledException);
}

void UM_constructor_stdstr(digest::UM_Digester& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h, uint64_t mod, uint64_t congruence){
	base_constructor_stdstr(dig, str, k, pos, minimized_h);
	CHECK(dig.get_mod() ==  mod);
	CHECK(dig.get_congruence() == congruence);
}

void base_constructor_cstr(digest::Digester& dig, const char* str, size_t len, unsigned k, size_t pos, unsigned minimized_h){
	CHECK(strcmp(str, dig.get_sequence()) == 0);
	CHECK(len == dig.get_len());
	CHECK(dig.get_k() == k);
	CHECK(dig.get_pos() == pos);
	CHECK(dig.get_minimized_h() == minimized_h);
	CHECK(dig.get_rolled() == false);
	CHECK_THROWS_AS(dig.get_chash(), digest::NotRolledException);
	CHECK_THROWS_AS(dig.get_fhash(), digest::NotRolledException);
	CHECK_THROWS_AS(dig.get_rhash(), digest::NotRolledException);
}

void UM_constructor_cstr(digest::UM_Digester& dig, const char* str, size_t len, unsigned k, size_t pos, unsigned minimized_h, uint64_t mod, uint64_t congruence){
	base_constructor_cstr(dig, str, len, k, pos, minimized_h);
	CHECK(dig.get_mod() ==  mod);
	CHECK(dig.get_congruence() == congruence);
}

void fhash_roll_one(digest::Digester& dig, std::string& str, unsigned k){
	nthash::NtHash tHash(str, 1, k, 0);
	uint64_t trueHash;
	uint64_t digHash;
	for(size_t i =0; i+k <= str.size(); i++){
		tHash.roll();
		trueHash = tHash.get_forward_hash();
		/*
		if(i == 0){
			trueHash= nthash::ntf64(str2.c_str(), k);
		}else{
			trueHash = nthash::ntf64(trueHash, k, str2[i-1], str2[i+k-1]);
		}
		*/
		dig.roll_one();
		digHash = dig.get_fhash();
		CHECK(dig.get_pos() == i);
		//std::cout << trueHash << " " << digHash << std::endl;
		//std::cout << i << " " << tHash.get_pos() << std::endl;
		CHECK(trueHash == digHash);
		CHECK(dig.get_rolled() == true);
	}

	CHECK_THROWS_AS(dig.roll_one(), std::out_of_range);
}

TEST_CASE("UM_Digester Testing"){
	setupStrings();
	/*
	for(unsigned i =0; i < strs.size(); i++){
		std::cout << strs[i] << std::endl;
	}
	*/
	SECTION("Testing Constructors"){
		unsigned k, minimized_h;
		uint64_t mod, congruence;
		size_t pos, len;
		std::string str;
		// string is length 1, k = 1
		str = "A";
		len = 1;
		k = ks[0];
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = i;
			mod = 2;
			congruence = 1;
			digest::UM_Digester* dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h);
			UM_constructor_stdstr(*dig, str, k, pos, minimized_h, mod, congruence);
			delete dig;

			dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h);
			UM_constructor_cstr(*dig, str.c_str(), len, k, pos, minimized_h, mod, congruence);
			delete dig;
		}
		// Using string in random.txt
		len = test_strs[4].size();
		k = ks[4];
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = i;
			mod = 1e9+7;
			congruence = 0;
			digest::UM_Digester* dig = new digest::UM_Digester(test_strs[4], k, mod, congruence, pos, minimized_h);
			UM_constructor_stdstr(*dig, test_strs[4], k, pos, minimized_h, mod, congruence);
			delete dig;

			dig = new digest::UM_Digester(test_strs[4].c_str(), len, k, mod, congruence, pos, minimized_h);
			UM_constructor_cstr(*dig, test_strs[4].c_str(), len, k, pos, minimized_h, mod, congruence);
			delete dig;
		}

		// pos = len-k
		str = "ACTGACTG";
		len = 8;
		k = ks[1];
		pos = 4;
		for(int i =0; i < 3; i++){
			minimized_h = i;
			mod = 1e9+7;
			congruence = 0;
			digest::UM_Digester* dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h);
			UM_constructor_stdstr(*dig, str, k, pos, minimized_h, mod, congruence);
			delete dig;

			dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h);
			UM_constructor_cstr(*dig, str.c_str(), len, k, pos, minimized_h, mod, congruence);
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
		digest::UM_Digester* dig;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);

		k = 2;
		// pos > seq.size()
		pos = 9;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);

		pos = 0;
		// pos + k > seq.size()
		pos = 7;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);

		pos = 0;
		// minimized_h > 2
		minimized_h = 3;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);

		minimized_h = 0;
		// mod >= congruence
		mod = 2;
		congruence = 2;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadModException);
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str.c_str(), len, k, mod, congruence, pos, minimized_h), digest::BadModException);

		mod = 1e9+7;
		congruence = 0;
	}

	SECTION("Testing roll_one"){
		
		for(int i =0; i < test_strs.size(); i++){
			for(int j =0; j < 5; j++){
				//std::cout << i << " " << j << std::endl;
				digest::UM_Digester* dig = new digest::UM_Digester(test_strs[i], ks[j], 1e9+7, 0, 0, 1);
				fhash_roll_one(*dig, test_strs[i], ks[j]);
			}
		}
		
		

		
	}
}

















/*
TEST_CASE("Testing: ntHash can roll()"){
	std::string str = "tcagcgaacgtaactg";
	std::cout << str << std::endl;
	nthash::NtHash test_hash(str, 1, 4);
	REQUIRE(test_hash.roll() == 1);
	//REQUIRE(test_hash->roll() == 0);
}

TEST_CASE("Quick Testing of UM_Digester"){
	std::string str = "TCAGTCTGAGACTCTGAGAGTCTGAGAGCTCTCTCGGGGTGTGTGCTA";
	std::cout << str << std::endl;
	
	digest::UM_Digester test_UM(str, 4, 2, 0, 0, 1);
	std::cout << test_UM.get_minimized_h() << std::endl;
	//test_UM.get_fhash();
	test_UM.roll_one();
	REQUIRE(test_UM.get_pos() == 0);
	test_UM.roll_next_minimizer();
	std::cout << test_UM.get_pos() << " " << test_UM.get_fhash() << " " << test_UM.get_chash() << std::endl;
	REQUIRE(test_UM.get_fhash() % 2 ==0);

}

TEST_CASE("Quick Testing of Replacing Strings"){
	std::string str = "TCAGTCTGAGACTCTGAGAGTCTGAGAGCTCTCTCGGGGTGTGTGCTA";
	std::cout << str << std::endl;
	
	digest::UM_Digester test_UM(str, 4, 2, 0, 0, 1);
	std::cout << test_UM.get_minimized_h() << std::endl;
	//test_UM.get_fhash();
	test_UM.roll_one();
	uint64_t og_chash = test_UM.get_chash();
	uint64_t og_fhash = test_UM.get_fhash();
	test_UM.roll_next_minimizer();
	unsigned og_pos = test_UM.get_pos();
	uint64_t og_min_fhash = test_UM.get_fhash();

	std::string str2 = "TCAGTCTGAGACTCTGAGAGTCTGAGAGCTCTCTCGGGGTGTGTGCTA";
	test_UM.new_seq(str2, 0);

	test_UM.roll_one();
	REQUIRE(test_UM.get_chash() == og_chash);
	REQUIRE(test_UM.get_fhash() == og_fhash);

	test_UM.roll_next_minimizer();

	REQUIRE(test_UM.get_pos() == og_pos);
	REQUIRE(test_UM.get_fhash() == og_min_fhash);



}

TEST_CASE("Quick Testing of Appending Strings"){
	std::string str = "TCAG";
	std::cout << str << std::endl;
	
	digest::UM_Digester test_UM(str, 4, 2, 0, 0, 1);
	std::cout << test_UM.get_minimized_h() << std::endl;
	//test_UM.get_fhash();
	test_UM.roll_one();
	
	REQUIRE(test_UM.roll_next_minimizer() == false);

	std::string str2 = "AAACCCTTTGGG";
	test_UM.append_seq(str2);

	test_UM.get_fhash();

	test_UM.roll_one();
	REQUIRE(test_UM.get_pos() == 1);

	test_UM.roll_next_minimizer();

	REQUIRE(test_UM.get_fhash() % 2 == 0);
	

}
*/