#include "digest/syncmer.hpp"

namespace digest{
	template <class T>
    void Syncmer<T>::roll_minimizer(unsigned amount, std::vector<uint32_t>& vec){
		amount += vec.size();

		while (this->st_size + 1 < this->large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				this->ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				this->ds.insert(this->get_pos(), this->fhash);
			}else{
				this->ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			this->st_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            Syncmer::fill_st(vec);
        }
    }

	template <class T>
    void Syncmer<T>::roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec){
		amount += vec.size();

		while (this->st_size + 1 < this->large_window and this->is_valid_hash) {
			if(this->get_minimized_h() == digest::MinimizedHashType::CANON){
				this->ds.insert(this->get_pos(), this->chash);
			}else if(this->get_minimized_h() == digest::MinimizedHashType::FORWARD){
				this->ds.insert(this->get_pos(), this->fhash);
			}else{
				this->ds.insert(this->get_pos(), this->rhash);
			}
			
			this->roll_one();
			this->st_size++;
		}

        while (this->is_valid_hash and vec.size() < amount){
            Syncmer::fill_st(vec);
        }
    }     

	template<class T>
    void Syncmer<T>::fill_st(std::vector<uint32_t>& vec){
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

	template<class T>
    void Syncmer<T>::fill_st(std::vector<std::pair<uint32_t, uint32_t>>& vec){
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
