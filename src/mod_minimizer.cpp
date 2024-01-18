#include "digest/mod_minimizer.hpp"
#include <cstdint>

namespace digest{
    void ModMin::roll_minimizer(unsigned amount, std::vector<uint32_t>& vec){
        if(!is_valid_hash) return;

		if(get_minimized_h() == digest::MinimizedHashType::CANON) {
			do {
                if((uint32_t)chash % mod == congruence){
                    vec.emplace_back(get_pos());
                }
			} while(roll_one() && vec.size() < amount);
			return;
		}

        if(get_minimized_h() == digest::MinimizedHashType::FORWARD) {
			do {
                if((uint32_t)fhash % mod == congruence){
                    vec.emplace_back(get_pos());
                }
			} while(roll_one() && vec.size() < amount);
			return;
		}

		// reverse
		do {
			if((uint32_t)rhash % mod == congruence){
				vec.emplace_back(get_pos());
			}
        } while(roll_one() && vec.size() < amount);
    }

	void ModMin::roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec){
        if(!is_valid_hash) return;

		if(get_minimized_h() == digest::MinimizedHashType::CANON) {
			do {
                if((uint32_t)chash % mod == congruence){
                    vec.emplace_back(get_pos(), chash);
                }
			} while(roll_one() && vec.size() < amount);
			return;
		}

        if(get_minimized_h() == digest::MinimizedHashType::FORWARD) {
			do {
                if((uint32_t)fhash % mod == congruence){
                    vec.emplace_back(get_pos(), fhash);
                }
			} while(roll_one() && vec.size() < amount);
			return;
		}

		// reverse
		do {
			if((uint32_t)rhash % mod == congruence){
				vec.emplace_back(get_pos(), rhash);
			}
        } while(roll_one() && vec.size() < amount);
    }
}
