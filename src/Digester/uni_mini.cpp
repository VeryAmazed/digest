#include "uni_mini.hpp"

namespace digest{
    bool UM_Digester::roll_next_minimizer(){
        while(end < len){
            roll_one();
            if(get_minimized_h() == 0){
                if(chash % mod == congruence) return true;
            }else if(get_minimized_h() == 1){
                if(fhash % mod == congruence) return true;
            }else{
                if(rhash % mod == congruence) return true;
            }
        }
        return false;
    }
}
