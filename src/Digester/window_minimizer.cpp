#include "window_minimizer.hpp"

namespace digest{
    void WindowMin::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        while(is_valid_hash){
            while(st_size < large_wind_kmer_am && is_valid_hash){
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
            if(st_size != large_wind_kmer_am){
                return;
            }
            st_size--;
            if(is_minimized){
                //std::cout << (*(st->segtree))[1].first << " " << (*(st->segtree))[1].second << " " << prev_mini.first << " " << prev_mini.second << std::endl;
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
}