#include "nthash.hpp"
#include <catch2/catch_test_macros.hpp>
#include "uni_mini.hpp"
#include <fstream>

std::vector<std::string> strs;
unsigned ks[] = {1, 4, 9, 89, 128};

void setupStrings(){
	std::string str;
	std::fstream fs;

	fs.open("../test_strings/A.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/N.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();
	
	fs.open("../test_strings/G.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/T.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/C.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/a_lowercase.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/coronavirus.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/random.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/random_lowercase.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/salmonella_enterica.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();

	fs.open("../test_strings/salmonella_lowercase.txt", std::fstream::in);
	fs >> str;
	strs.push_back(str);
	fs.close();
}

void base_constructor_stdstr(digest::Digester& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h){
	CHECK(strcmp(str.c_str(), dig.get_sequence()) == 0);
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

TEST_CASE("UM_Digester Testing"){
	setupStrings();
	/*
	for(unsigned i =0; i < strs.size(); i++){
		std::cout << strs[i] << std::endl;
	}
	*/
	SECTION("Testing Constructors"){
		// k = 1, str_len = 1
		// normal
		// start at len-k
		// all the errors
		unsigned k, minimized_h;
		uint64_t mod, congruence;
		size_t pos;
		std::string str;
		// string is length 1, k = 1
		str = "A";
		k = 1;
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = i;
			mod = 2;
			congruence = 1;
			
		}
		// Using all A's
		for(int i =0; i < 3; i++){

		}
		// Using all a's
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