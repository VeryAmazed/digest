#include "digest/window_minimizer.hpp"

namespace digest{
	template <int32_t large_window>
    void WindowMin<large_window>::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        while(is_valid_hash){
            // -------------------
            fill_st();
            // -------------------
            if((st_size != large_window) || (vec.size() >= amount)){
                return;
            }
            st_size--;
            // -------------------
            check(vec);
        }
        
    }

	template <int32_t large_window>
    void WindowMin<large_window>::fill_st(){
        while((st_size < large_window) && is_valid_hash){
            if(get_minimized_h() == digest::MinimizedHashType::CANON){
                st.set(chash, get_pos());
            }else if(get_minimized_h() == digest::MinimizedHashType::FORWARD){
                st.set(fhash, get_pos());
            }else{
                st.set(rhash, get_pos());
            }
            st_size++;
            
            roll_one();
        }
    }

	template <int32_t large_window>
    void WindowMin<large_window>::check(std::vector<size_t>& vec){
        if(is_minimized){
            if(st.min() != prev_mini){
                prev_mini =	st.min();
                vec.emplace_back(prev_mini);
            }
        }else{
            is_minimized = true;
            prev_mini = st.min();
			vec.emplace_back(prev_mini);
        }
    }
}
