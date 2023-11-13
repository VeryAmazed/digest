#include "digest/syncmer.hpp"

namespace digest{
    void Syncmer::check1(std::vector<size_t>& vec){
        size_t right = (st_index + large_wind_kmer_am - 1) % large_wind_kmer_am;
        if(st.segtree[1].first == st.segtree[st_index + large_wind_kmer_am].first || st.segtree[1].first == st.segtree[right + large_wind_kmer_am].first){
            vec.push_back(st.segtree[st_index + large_wind_kmer_am].second);
        }
    }

    void Syncmer::check2(std::vector<std::pair<size_t, size_t>>& vec){
        size_t right = (st_index + large_wind_kmer_am - 1) % large_wind_kmer_am;
        if(st.segtree[1].first == st.segtree[st_index + large_wind_kmer_am].first || st.segtree[1].first == st.segtree[right + large_wind_kmer_am].first){
            vec.push_back(std::make_pair(st.segtree[st_index + large_wind_kmer_am].second,  st.segtree[right + large_wind_kmer_am].second));
            
        }
    }

    void Syncmer::roll_minimizer(unsigned amount, std::vector<size_t>& vec){
        while(is_valid_hash){
            // -------------------
            fill_st();
            // -------------------
            if((st_size != large_wind_kmer_am) || (vec.size() >= amount)){
                return;
            }
            st_size--;
            // -------------------
            check1(vec);
        }
    }

    void Syncmer::roll_minimizer(unsigned amount, std::vector<std::pair<size_t, size_t>>& vec){
        while(is_valid_hash){
            // -------------------
            fill_st();
            // -------------------
            if((st_size != large_wind_kmer_am) || (vec.size() >= amount)){
                return;
            }
            st_size--;
            // -------------------
            check2(vec);
        }
    }     
}
