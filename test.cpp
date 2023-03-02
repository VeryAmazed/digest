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
	nthash::NtHash* test_hash = new nthash::NtHash(str, 1, 4);
	REQUIRE(test_hash->roll() == 1);
	//REQUIRE(test_hash->roll() == 0);
}