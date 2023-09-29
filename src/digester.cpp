#include "digest/digester.hpp"

namespace digest{
    

    void Digester::append_seq(const char* seq, size_t len){
        if(end < this->len){
            throw NotRolledTillEndException();
        }
        offset += this->len;
        size_t ind = this->len-1;
        
        /*
            this is for the case where we call append_seq after having previously called append_seq and not having gotten through the deque
            In such a case, since append_seq initializes a hash, we need to get rid of the first character in the deque since if we just initialized the hash
            without doing this, it would be identical the the current hash held by the object
            
            However, there is also the case that a hash was never previously initialized, such as when the length of the string used in the previous append_seq call, 
            plus the the amount of ACTG characters after the last non-ACTG character in the original string summed to be less than k
            In this case, it would not be correct to remove the first character in the deque
        */
        if((start != end || c_outs->size() == k) && c_outs->size() > 0){
            c_outs->pop_front();
        }

        // the following copies in characters from the end of the old sequence into the deque
        std::vector<char> temp_vec;
        while(temp_vec.size() + c_outs->size()< k-1 && ind >= start){
            if(!is_ACTG(this->seq[ind])) break;
            
            temp_vec.push_back(this->seq[ind]);
            if(ind == 0) break;

            ind--;
        }
        for(std::vector<char>::reverse_iterator rit = temp_vec.rbegin(); rit != temp_vec.rend(); rit++){
            c_outs->push_back(*rit);
        }

        // the following copies in characters from the front of the new sequence if there weren't enough non-ACTG characters at the end of the old sequence
        ind = 0;
        start = 0;
        end = 0;
        while(c_outs->size() < k && ind < len){
            if(!is_ACTG(seq[ind])){
                start = ind+1;
                end = start + k;
                this->seq = seq;
                this->len = len;
                c_outs->clear();
                init_hash();
                break;
            }
            c_outs->push_back(seq[ind]);
            ind++;
            start++;
            end++;
        }

        // the following initializes a hash if we managed to fill the deque
        if(c_outs->size() == k){
            std::string temp(c_outs->begin(), c_outs->end());
            unsigned locn_useless;
            // nthash::ntc64(temp.c_str(), k, fhash, rhash, chash, locn_useless);
            is_valid_hash = true;

        }
        this->seq = seq;
        this->len = len;
    }

    void Digester::append_seq(const std::string& seq){
        append_seq(seq.c_str(), seq.size());
    }

    bool Digester::init_hash(){
        c_outs->clear();
        unsigned locn_useless;
        while(end-1 < len){
            bool works = true;
            for(size_t i = start; i < end; i++){
                if(!is_ACTG(seq[i])){
                    start = i+1;
                    end = start + k;
                    works = false;
                    break;
                }
            }
            if(!works){
                continue;
            }
            // nthash::ntc64(seq + start, k, fhash, rhash, chash, locn_useless);
            is_valid_hash = true;
            return true;
        }
        is_valid_hash = false;
        return false;
    }

    bool Digester::roll_one(){
        if(!is_valid_hash){
            return false;
        }
        if(end >= len){
            is_valid_hash = false;
            return false;
        }
        if(c_outs->size() > 0){
            if(is_ACTG(seq[end])){
                // fhash = nthash::ntf64(fhash, k, c_outs->front(), seq[end]);
                // rhash = nthash::ntr64(rhash, k, c_outs->front(), seq[end]);
                c_outs->pop_front();
                end++;
                // chash = nthash::canonical(fhash,rhash);
                return true;
            }else{
                // c_outs will contain at most k-1 characters, so if we jump to end + 1, we won't consider anything else in deque so we should clear it
                c_outs->clear();
                start = end+1;
                end = start + k;
                return init_hash();
            }
        }else{
            if(is_ACTG(seq[end])){
                // fhash = nthash::ntf64(fhash, k, seq[start], seq[end]);
                // rhash = nthash::ntr64(rhash, k, seq[start], seq[end]);
                start++;
                end++;
                // chash = nthash::canonical(fhash,rhash);
                return true;
            }else{
                start = end+1;
                end = start + k;
                return init_hash();
            }
        }
        
    }
} 
