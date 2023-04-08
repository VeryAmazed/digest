#include "nthash.hpp"
#include <catch2/catch_test_macros.hpp>
#include "uni_mini.hpp"
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

void base_constructor_stdstr(digest::Digester& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h){
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

void UM_constructor_stdstr(digest::UM_Digester& dig, std::string& str, unsigned k, size_t pos, unsigned minimized_h, uint64_t mod, uint64_t congruence){
	base_constructor_stdstr(dig, str, k, pos, minimized_h);
	CHECK(dig.get_mod() ==  mod);
	CHECK(dig.get_congruence() == congruence);
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

void um_roll_minimizer(digest::UM_Digester& dig, std::string& str, unsigned k, unsigned minimized_h, uint64_t prime){
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
	std::vector<size_t> dig_positions = dig.roll_minimizer(400);
	REQUIRE(positions.size() == dig_positions.size());
	for(size_t i = 0; i < positions.size(); i++){
		CHECK(dig_positions[i] == positions[i]);
	}
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
		}

		// string is length 1, k = 4
		str = "A";
		len = 1;
		k = ks[1];
		pos = 0;
		for(int i =0; i < 3; i++){
			minimized_h = i;
			mod = 2;
			congruence = 1;
			digest::UM_Digester* dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h);
			UM_constructor_stdstr(*dig, str, k, pos, minimized_h, mod, congruence);
			delete dig;
		}	
		
		for(int i =0; i < test_strs.size(); i++){
			len = test_strs[i].size();
			for(int j =0; j < 8; j++){
				k = ks[j];
				for(int l =0; l < 16; l++){
					pos = l;
					for(int p =0; p < 3; p++){
						minimized_h = p;
						mod = 1e9+7;
						congruence = 0;
						digest::UM_Digester* dig = new digest::UM_Digester(test_strs[i], k, mod, congruence, pos, minimized_h);
						UM_constructor_stdstr(*dig, test_strs[i], k, pos, minimized_h, mod, congruence);
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
		mod = 1e9+7;
		congruence = 0;
		// k = 0
		k = 0;
		digest::UM_Digester* dig;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		k = 2;

		// pos >= seq.size()
		pos = 8;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		pos = 0;

		// minimized_h > 2
		minimized_h = 3;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadConstructionException);
		minimized_h = 0;

		// mod >= congruence
		mod = 2;
		congruence = 2;
		CHECK_THROWS_AS(dig = new digest::UM_Digester(str, k, mod, congruence, pos, minimized_h), digest::BadModException);
		mod = 1e9+7;
		congruence = 0;
	}

	
	SECTION("Testing roll_one"){
		for(int i =0; i < 7; i++){
			for(int j =0; j < 8; j++){
				digest::UM_Digester* dig = new digest::UM_Digester(test_strs[i], ks[j], 1e9+7, 0, 0, 1);
				roll_one(*dig, test_strs[i], ks[j]);
				delete dig;
			}
		}
	}
	
	SECTION("Testing roll_minimizer(). The one that takes no parameters"){
		uint64_t prime = 17;
		for(int i =0; i < 7; i += 2){
			for(int j =0; j < 8; j++){
				for(int l = 0; l < 3; l++){
					digest::UM_Digester* dig = new digest::UM_Digester(test_strs[i], ks[j], prime, 0, 0, l);
					um_roll_minimizer(*dig, test_strs[i], ks[j], l, prime);
					delete dig;
				}
				
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