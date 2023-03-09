#include "digester.hpp"

namespace digest{
    void Digester::append_seq(const std::string& seq){
        if(!rolled){
            throw NotRolledException();
        }
        size_t ind = end;
        while(c_outs->size() < k){
            c_outs->push_front(this->seq[ind]);
            ind--;
        }
        start = 0;
        end =0;
        this->seq = seq.data();
        this->len = seq.size();
    }

    void Digester::append_seq(const char* seq, size_t len){
        if(!rolled){
            throw NotRolledException();
        }
        size_t ind = end;
        while(c_outs->size() < k){
            c_outs->push_front(this->seq[ind]);
            ind--;
        }
        start = 0;
        end =0;
        this->seq = seq;
        this->len = len;
    }

    void Digester::roll_one(){
        
        if(!rolled){
            if(end-1 == len){
                throw std::out_of_range("End of sequence.");
            }
            fhash = nthash::ntf64(seq, k);
            rhash = nthash::ntr64(seq, k);
            rolled = true;
        }else{
            if(end == len){
                throw std::out_of_range("End of sequence.");
            }
            if(c_outs->size() > 0){
                fhash = nthash::ntf64(fhash, k, seq[c_outs->front()], seq[end]);
                rhash = nthash::ntr64(rhash, k, seq[c_outs->front()], seq[end]);
                c_outs->pop_front();
                /*
                if(c_outs->size() == 0){
                    start = 0;
                }
                */
            }else{
                fhash = nthash::ntf64(fhash, k, seq[start], seq[end]);
                rhash = nthash::ntr64(rhash, k, seq[start], seq[end]);
                start++;
            }
            end++;
            pos++;
        }
        chash = nthash::canonical(fhash,rhash);
        
    }

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
} 
