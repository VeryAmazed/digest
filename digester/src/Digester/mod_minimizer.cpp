#include "mod_minimizer.hpp"

namespace digest{
    void ModMin::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        if(!is_valid_hash) return;
        do{
            if(get_minimized_h() == 0){
                if(chash % mod == congruence){
                    vec.push_back(get_pos());
                } 
            }else if(get_minimized_h() == 1){
                if(fhash % mod == congruence){
                    vec.push_back(get_pos());
                }
            }else{
                if(rhash % mod == congruence){
                    vec.push_back(get_pos());
                }
            }
        }while(roll_one() && vec.size() < amount);

    }
}
