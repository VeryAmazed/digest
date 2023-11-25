#include "digest/mod_minimizer.hpp"

namespace digest{
    void ModMin::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        if(!is_valid_hash) return;

		if(get_minimized_h() == digest::MinimizedHashType::CANON) {
			do {
                if(chash % mod == congruence){
                    vec.push_back(get_pos());
                }
			} while(roll_one() && vec.size() < amount);
			return;
		}

        if(get_minimized_h() == digest::MinimizedHashType::FORWARD) {
			do {
                if(fhash % mod == congruence){
                    vec.push_back(get_pos());
                }
			} while(roll_one() && vec.size() < amount);
			return;
		}

		do {
			if(rhash % mod == congruence){
				vec.push_back(get_pos());
			}
        } while(roll_one() && vec.size() < amount);
    }
}
