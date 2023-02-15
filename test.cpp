#include "nthash.hpp"

int main(){
	std::string str = "TGACGCATTTGGCGAGATGATG";
	nthash::NtHash* test_hash = new nthash::NtHash(str, 1, 4);
	std::cout << test_hash->roll() << std::endl;
	
	return 0;
}