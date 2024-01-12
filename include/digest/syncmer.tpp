#include "digest/syncmer.hpp"

namespace digest{

//	template <class T>
//    void Syncmer<T>::check2(std::vector<std::pair<size_t, size_t>>& vec){
//        //size_t right = (st.i + this->large_window - 1) % this->large_window;
//        size_t right = (this->st.i + this->large_window - 1);
//        if(right >= this->large_window){
//            right -= this->large_window;
//        }
//        if(this->st.get_hash(1) == this->st.get_hash(this->st.i + this->large_window) || this->st.get_hash(1) == this->st.get_hash(right + this->large_window)){
//            vec.push_back(std::make_pair(this->st.get_index(this->st.i + this->large_window),  this->st.get_index(right + this->large_window)));
//        }
//    }

	template <class T>
    void Syncmer<T>::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
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
            this->fill_st(vec);
        }
    }

//	template <class T>
//    void Syncmer<T>::roll_minimizer(unsigned amount, std::vector<std::pair<size_t, size_t>>& vec){
//        while(this->is_valid_hash){
//            // -------------------
//            this->fill_st(vec);
//            // -------------------
//            if((this->st_size != this->large_window) || (vec.size() >= amount)){
//                return;
//            }
//            this->st_size--;
//            // -------------------
//            check2(vec);
//        }
//    }     
}
