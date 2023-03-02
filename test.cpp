#include "nthash.hpp"
#include <catch2/catch_test_macros.hpp>
#include "uni_mini.hpp"
/*
int main(){
	std::string str = "TGACGCATTTGGCGAGATGATG";
	nthash::NtHash* test_hash = new nthash::NtHash(str, 1, 4);
	std::cout << test_hash->roll() << std::endl;
	
	return 0;
}
*/

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