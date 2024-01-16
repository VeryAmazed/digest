#include "digest/window_minimizer.hpp"

namespace digest{
	template<class T>
    void WindowMin<T>::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
		amount += vec.size();

		while (st_size + 1 < large_window and is_valid_hash) {
			if(get_minimized_h() == digest::MinimizedHashType::CANON){
				ds.insert(get_pos(), chash);
			}else if(get_minimized_h() == digest::MinimizedHashType::FORWARD){
				ds.insert(get_pos(), fhash);
			}else{
				ds.insert(get_pos(), rhash);
			}
			
			roll_one();
			st_size++;
		}

        while (is_valid_hash and vec.size() < amount){
            fill_st(vec);
        }
    }

	template<class T>
    void WindowMin<T>::fill_st(std::vector<size_t>& vec){
		if(get_minimized_h() == digest::MinimizedHashType::CANON){
			check(vec, ds.insert(get_pos(), chash));
		}else if(get_minimized_h() == digest::MinimizedHashType::FORWARD){
			check(vec, ds.insert(get_pos(), fhash));
		}else{
			check(vec, ds.insert(get_pos(), rhash));
		}
		
		roll_one();
    }

	template<class T>
    void WindowMin<T>::check(std::vector<size_t>& vec, uint32_t ind){
        if(is_minimized){
            if(ind != prev_mini){
                prev_mini =	ind;
                vec.emplace_back(prev_mini);
            }
        }else{
            is_minimized = true;
            prev_mini = ind;
			vec.emplace_back(prev_mini);
        }
    }
}
