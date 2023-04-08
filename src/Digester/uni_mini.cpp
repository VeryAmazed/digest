#include "uni_mini.hpp"

namespace digest{
    std::vector<size_t> UM_Digester::roll_minimizer(unsigned amount){
        std::vector<size_t> vec;
        if(!is_valid_hash) return vec;
        do{
            if(get_minimized_h() == 0){
                if(chash % mod == congruence){
                    vec.push_back(pos);
                } 
            }else if(get_minimized_h() == 1){
                if(fhash % mod == congruence){
                    vec.push_back(pos);
                }
            }else{
                if(rhash % mod == congruence){
                    vec.push_back(pos);
                }
            }
        }while(roll_one() && vec.size() < amount);
        return vec;
    }
}
