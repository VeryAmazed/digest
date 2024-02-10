#include "digest/window_minimizer.hpp"

namespace digest{
	template<BadCharPolicy P, class T>
    void WindowMin<P, T>::roll_minimizer(unsigned amount, std::vector<uint32_t>& vec){
		amount += vec.size();

		while (ds_size + 1 < large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				ds.insert(this->get_pos(), this->fhash);
			}else{
				ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			ds_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            roll_ds_wind(vec);
        }
    }

	template<BadCharPolicy P, class T>
    void WindowMin<P, T>::roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec){
		amount += vec.size();

		while (ds_size + 1 < large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				ds.insert(this->get_pos(), this->fhash);
			}else{
				ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			ds_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            roll_ds_wind(vec);
        }
    }

	template<BadCharPolicy P, class T>
    void WindowMin<P, T>::roll_ds_wind(std::vector<uint32_t>& vec){
		if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
			ds.insert(this->get_pos(), this->chash);
		}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
			ds.insert(this->get_pos(), this->fhash);
		}else{
			ds.insert(this->get_pos(), this->rhash);
		}
		check(vec);
		
		this->roll_one();
    }

	template<BadCharPolicy P, class T>
    void WindowMin<P, T>::roll_ds_wind(std::vector<std::pair<uint32_t, uint32_t>>& vec){
		if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
			ds.insert(this->get_pos(), this->chash);
		}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
			ds.insert(this->get_pos(), this->fhash);
		}else{
			ds.insert(this->get_pos(), this->rhash);
		}
		check(vec);
		
		this->roll_one();
    }

	template<BadCharPolicy P, class T>
    void WindowMin<P, T>::check(std::vector<uint32_t>& vec){
        if(is_minimized){
            if(ds.min() != prev_mini){
                prev_mini =	ds.min();
                vec.emplace_back(prev_mini);
            }
        }else{
            is_minimized = true;
            prev_mini = ds.min();
			vec.emplace_back(prev_mini);
        }
    }

	template<BadCharPolicy P, class T>
    void WindowMin<P, T>::check(std::vector<std::pair<uint32_t, uint32_t>>& vec){
        if(is_minimized){
            if(ds.min() != prev_mini){
                prev_mini =	ds.min();
                vec.emplace_back(prev_mini, ds.min_hash());
            }
        }else{
            is_minimized = true;
            prev_mini = ds.min();
			vec.emplace_back(prev_mini, ds.min_hash());
        }
    }
}
