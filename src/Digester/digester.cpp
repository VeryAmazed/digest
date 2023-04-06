#include "digester.hpp"

namespace digest{
    

    void Digester::append_seq(const char* seq, size_t len){
        if(end >= this->len){
            throw NotRolledTillEndException();
        }
        size_t ind = end-1;
        if(is_valid_hash){
            while(c_outs->size() <= k && ind >= 0){
                c_outs->push_front(this->seq[ind]);
                ind--;
            }
            start = 0;
            end =0;
        }else{
            while(c_outs->size() <= k && ind >= 0){
                if(!is_ACTG(this->seq[ind])){
                    break;
                }
                c_outs->push_front(this->seq[ind]);
                ind--;
            }
            ind = 0;
            while(c_outs->size() <= k && ind < len){
                if(!is_ACTG(this->seq[ind])){
                    pos += c_outs->size() + 2;
                    start = ind+1;
                    end = start + k;
                    init_hash();
                    c_outs->clear();
                    break;
                }
                c_outs->push_back(seq[ind]);
                ind++;
                end++;
            }
            if(c_outs->size() == k){
                std::string temp(c_outs->begin(), c_outs->end());
                fhash = nthash::ntf64(temp.c_str(), k);
                rhash = nthash::ntr64(temp.c_str(), k);
            }
        }
        this->seq = seq;
        this->len = len;
    }

    void Digester::append_seq(const std::string& seq){
        append_seq(seq.c_str(), seq.size());
    }

    bool Digester::init_hash(){
        unsigned locn_useless;
        while(end-1 < len){
            bool works = true;
            for(size_t i = start; i < end; i++){
                if(!is_ACTG(seq[i])){
                    pos += (i+1) - start;
                    start = i+1;
                    end = start + k;
                    works = false;
                    break;
                }
            }
            if(!works){
                continue;
            }
            /*
            fhash = nthash::ntf64(seq + start, k);
            rhash = nthash::ntr64(seq + start, k);
            chash = nthash::canonical(fhash, rhash);
            */
            nthash::ntc64(seq + start, k, fhash, rhash, chash, locn_useless);
            is_valid_hash = true;
            return true;
        }
        is_valid_hash = false;
        return false;
    }

    bool Digester::roll_one(){
        if(!is_valid_hash){
            return false;
        }else{
            if(end >= len){
                is_valid_hash = false;
                return false;
            }
            if(c_outs->size() > 0){
                if(is_ACTG(seq[end])){
                    fhash = nthash::ntf64(fhash, k, seq[c_outs->front()], seq[end]);
                    rhash = nthash::ntr64(rhash, k, seq[c_outs->front()], seq[end]);
                    c_outs->pop_front();
                    pos++;
                    end++;
                    return true;
                }else{
                    pos += k+1;
                    c_outs->clear();
                    start = end+1;
                    end = start + k;
                    return init_hash();
                }
            }else{
                if(is_ACTG(seq[end])){
                    fhash = nthash::ntf64(fhash, k, seq[start], seq[end]);
                    rhash = nthash::ntr64(rhash, k, seq[start], seq[end]);
                    start++;
                    pos++;
                    end++;
                    return true;
                }else{
                    pos += k+1;
                    start = end+1;
                    end = start + k;
                    return init_hash();
                }
                
            }
            chash = nthash::canonical(fhash,rhash);
        }
    }
    /*
    std::string Digester::get_string(){
        std::string str;
        if(c_outs->size() > 0){
            std::deque<char>::iterator it = c_outs->begin();
            while(it != c_outs->end()){
                str.push_back(*it);
            }
        }
        for(int i = start; i < end; i++){
            str.push_back(seq[i]);
        }
        return str;
    }
    */
} 
