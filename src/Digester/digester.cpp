#include "digester.hpp"

namespace digest{
    void Digester::append_seq(const std::string& seq){
        if(!rolled){
            throw NotRolledException();
        }
        size_t ind = end;
        while(c_outs.size() < k){
            c_outs.push_front(this->seq[ind]);
            ind--;
        }

        end =0;
        this->seq = seq.data();
        this->len = seq.size();
    }

    void Digester::append_seq(const char* seq, size_t len){
        if(!rolled){
            throw NotRolledException();
        }
        size_t ind = end;
        while(c_outs.size() < k){
            c_outs.push_front(this->seq[ind]);
            ind--;
        }

        end =0;
        this->seq = seq;
        this->len = len;
    }

    void Digester::roll_one(){
        if(end == len){
            throw std::out_of_range("End of sequence.");
        }
        if(!rolled){
            fhash = nthash::ntf64(seq, k);
            rhash = nthash::ntr64(seq, k);
            rolled = true;
        }else{
            if(c_outs.size() > 0){
                fhash = nthash::ntf64(fhash, k, seq[c_outs.front()], seq[end]);
                rhash = nthash::ntr64(rhash, k, seq[c_outs.front()], seq[end]);
                c_outs.pop_front();
                if(c_outs.size() == 0){
                    start = 0;
                }
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
} 
