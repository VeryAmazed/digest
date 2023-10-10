#include "digest/window_minimizer.hpp"

namespace digest{
    void WindowMin::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        while(is_valid_hash){
            // -------------------
            fill_st();
            // -------------------
            if((st_size != large_wind_kmer_am) || (vec.size() >= amount)){
                return;
            }
            st_size--;
            // -------------------
            check(vec);
        }
        
    }

    void WindowMin::fill_st(){
        while((st_size < large_wind_kmer_am) && is_valid_hash){
            if(get_minimized_h() == 0){
                st->set(st_index, std::make_pair(chash, get_pos()));
            }else if(get_minimized_h() == 1){
                st->set(st_index, std::make_pair(fhash, get_pos()));
            }else{
                st->set(st_index, std::make_pair(rhash, get_pos()));
            }
            st_size++;
            st_index = (st_index + 1) % large_wind_kmer_am;
            roll_one();
        }
    }

    void WindowMin::check(std::vector<size_t>& vec){
        if(is_minimized){
            if((*(st->segtree))[1] != prev_mini){
                prev_mini = (*(st->segtree))[1];
                vec.push_back((*(st->segtree))[1].second);
            }
        }else{
            is_minimized = true;
            prev_mini = (*(st->segtree))[1];
            vec.push_back((*(st->segtree))[1].second);
        }
    }
}
