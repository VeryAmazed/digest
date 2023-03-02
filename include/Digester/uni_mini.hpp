#ifndef UNI_MINI_HPP
#define UNI_MINI_HPP

#include "digester.hpp"

namespace digest{

class UM_Digester : public Digester{
    public:
        UM_Digester(const std::string& seq, unsigned k, uint64_t mod, uint64_t congruence, size_t pos = 0, bool canonicalized = true) :  Digester(seq, k, pos, canonicalized){
           this->mod = mod;
           this->congruence = congruence;
        }

        std::pair<bool, size_t> roll() override;

        std::string get_min_string() override;

        uint64_t get_min_hash() = 0;

        void change_seq(const std::string& new_seq, size_t new_pos) override;
    
    private:
        uint64_t mod;
        uint64_t congruence;

};

}

#endif