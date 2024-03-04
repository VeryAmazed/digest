#include "digest/mod_minimizer.hpp"
#include <cstdint>

namespace digest{

	template <BadCharPolicy P>
    void ModMin<P>::roll_minimizer(unsigned amount, std::vector<uint32_t>& vec){
        if(!this->is_valid_hash) return;

		if(this->get_minimized_h() == digest::MinimizedHashType::CANON) {
			do {
                if((uint32_t)this->chash % mod == congruence){
                    vec.emplace_back(this->get_pos());
                }
			} while(this->roll_one() && vec.size() < amount);
			return;
		}

        if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD) {
			do {
                if((uint32_t)this->fhash % mod == congruence){
                    vec.emplace_back(this->get_pos());
                }
			} while(this->roll_one() && vec.size() < amount);
			return;
		}

		// reverse
		do {
			if((uint32_t)this->rhash % mod == congruence){
				vec.emplace_back(this->get_pos());
			}
        } while(this->roll_one() && vec.size() < amount);
    }

	template <BadCharPolicy P>
	void ModMin<P>::roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec){
        if(!this->is_valid_hash) return;

		if(this->get_minimized_h() == digest::MinimizedHashType::CANON) {
			do {
                if((uint32_t)this->chash % mod == congruence){
                    vec.emplace_back(this->get_pos(), this->chash);
                }
			} while(this->roll_one() && vec.size() < amount);
			return;
		}

        if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD) {
			do {
                if((uint32_t)this->fhash % mod == congruence){
                    vec.emplace_back(this->get_pos(), this->fhash);
                }
			} while(this->roll_one() && vec.size() < amount);
			return;
		}

		// reverse
		do {
			if((uint32_t)this->rhash % mod == congruence){
				vec.emplace_back(this->get_pos(), this->rhash);
			}
        } while(this->roll_one() && vec.size() < amount);
    }
}
