#include "digester.hpp"

namespace digest{
    

    void Digester::append_seq(const char* seq, size_t len){
        if(end < this->len){
            throw NotRolledTillEndException();
        }
        offset += this->len;
        size_t ind = this->len-1;
        //std::cout << ind <<std::endl;
        
        while(c_outs->size() < k-1 && ind >= 0){
            if(!is_ACTG(this->seq[ind])){
                break;
            }
            c_outs->push_front(this->seq[ind]);
            //std::cout << c_outs->front() <<std::endl;
            ind--;
        }
        ind = 0;
        
        start = 0;
        end = 0;
        while(c_outs->size() < k && ind < len){
            if(!is_ACTG(seq[ind])){
                start = ind+1;
                end = start + k;
                this->seq = seq;
                this->len = len;
                init_hash();
                c_outs->clear();
                break;
            }
            c_outs->push_back(seq[ind]);
            //std::cout << c_outs->back() <<std::endl;
            ind++;
            start++;
            end++;
        }
        if(c_outs->size() == k){
            std::string temp(c_outs->begin(), c_outs->end());
            unsigned locn_useless;
            nthash::ntc64(temp.c_str(), k, fhash, rhash, chash, locn_useless);
            is_valid_hash = true;

        }
        
        /*
        for(std::deque<char>::iterator it = c_outs->begin(); it != c_outs->end(); it++){
            std::cout << *it << std::endl;
        }
        */
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
                    start = i+1;
                    end = start + k;
                    works = false;
                    break;
                }
            }
            if(!works){
                continue;
            }
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
                    fhash = nthash::ntf64(fhash, k, c_outs->front(), seq[end]);
                    rhash = nthash::ntr64(rhash, k, c_outs->front(), seq[end]);
                    c_outs->pop_front();
                    end++;
                    chash = nthash::canonical(fhash,rhash);
                    return true;
                }else{
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
                    end++;
                    chash = nthash::canonical(fhash,rhash);
                    return true;
                }else{
                    start = end+1;
                    end = start + k;
                    return init_hash();
                }
                
            }
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
