#include "digest/syncmer.hpp"

namespace digest{

	template <int32_t large_window>
    void Syncmer<large_window>::check1(std::vector<size_t>& vec){
        //size_t right = (st.i + large_window - 1) % large_window;
        size_t right = (this->st.i + large_window - 1);
        if(right >= large_window){
            right -= large_window;
        }
        if(this->st.get_hash(1) == this->st.get_hash(this->st.i + large_window) || this->st.get_hash(1) == this->st.get_hash(right + large_window)){
            vec.push_back(this->st.get_index(this->st.i + large_window));
        }
    }

	template <int32_t large_window>
    void Syncmer<large_window>::check2(std::vector<std::pair<size_t, size_t>>& vec){
        //size_t right = (st.i + large_window - 1) % large_window;
        size_t right = (this->st.i + large_window - 1);
        if(right >= large_window){
            right -= large_window;
        }
        if(this->st.get_hash(1) == this->st.get_hash(this->st.i + large_window) || this->st.get_hash(1) == this->st.get_hash(right + large_window)){
            vec.push_back(std::make_pair(this->st.get_index(this->st.i + large_window),  this->st.get_index(right + large_window)));
        }
    }

	template <int32_t large_window>
    void Syncmer<large_window>::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        while(this->is_valid_hash){
            // -------------------
            this->fill_st();
            // -------------------
            if((this->st_size != large_window) || (vec.size() >= amount)){
                return;
            }
            this->st_size--;
            // -------------------
            check1(vec);
        }
    }

	template <int32_t large_window>
    void Syncmer<large_window>::roll_minimizer(unsigned amount, std::vector<std::pair<size_t, size_t>>& vec){
        while(this->is_valid_hash){
            // -------------------
            this->fill_st();
            // -------------------
            if((this->st_size != large_window) || (vec.size() >= amount)){
                return;
            }
            this->st_size--;
            // -------------------
            check2(vec);
        }
    }     
}
