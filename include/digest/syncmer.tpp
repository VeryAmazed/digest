#include "digest/syncmer.hpp"

namespace digest{
	template <BadCharPolicy P, class T>
    void Syncmer<P, T>::roll_minimizer(unsigned amount, std::vector<uint32_t>& vec){
		amount += vec.size();

		while (this->ds_size + 1 < this->large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				this->ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				this->ds.insert(this->get_pos(), this->fhash);
			}else{
				this->ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			this->ds_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            Syncmer::roll_ds_sync(vec);
        }
    }

	template <BadCharPolicy P, class T>
    void Syncmer<P, T>::roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec){
		amount += vec.size();

		while (this->ds_size + 1 < this->large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				this->ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				this->ds.insert(this->get_pos(), this->fhash);
			}else{
				this->ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			this->ds_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            Syncmer::roll_ds_sync(vec);
        }
    }     

	template<BadCharPolicy P, class T>
    void Syncmer<P, T>::roll_ds_sync(std::vector<uint32_t>& vec){
		if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
			this->ds.insert(this->get_pos(), this->chash);
		}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
			this->ds.insert(this->get_pos(), this->fhash);
		}else{
			this->ds.insert(this->get_pos(), this->rhash);
		}
		this->ds.min_syncmer(vec);
		
		this->roll_one();
    }

	template<BadCharPolicy P, class T>
    void Syncmer<P, T>::roll_ds_sync(std::vector<std::pair<uint32_t, uint32_t>>& vec){
		if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
			this->ds.insert(this->get_pos(), this->chash);
		}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
			this->ds.insert(this->get_pos(), this->fhash);
		}else{
			this->ds.insert(this->get_pos(), this->rhash);
		}
		this->ds.min_syncmer(vec);
		
		this->roll_one();
    }
}