#include "uni_mini.hpp"

namespace digest{
    void UM_Digester::roll_next_minimizer(){
        while(end < len){
            roll_one();
            if(get_minimized_h() == 0){
                if(chash % mod == congruence) break;
            }else if(get_minimized_h() == 1){
                if(fhash % mod == congruence) break;
            }else{
                if(rhash % mod == congruence) break;
            }
        }
    }
}
