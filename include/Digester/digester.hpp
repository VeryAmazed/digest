#ifndef DIGESTER_HPP
#define DIGESTER_HPP

#include "../ntHash/nthash.hpp"

namespace digest{
typedef long long ll;

#define pb push_back 
#define f first
#define s second
#define mp make_pair
#define pll pair<ll, ll>
#define pii pair<int, int>

class Digester{
    public:
        
        Digester(const std::string& seq, unsigned k, size_t pos = 0){
            hasher = new nthash::NtHash(seq, 1, k, pos);
        }

        virtual bool roll() = 0;

        virtual uint64_t get_canonical_m() = 0;

        virtual uint64_t get_forward_m() = 0;

        virtual uint64_t get_reverse_m() = 0;

        virtual void change_seq(const std::string& new_seq, size_t new_pos = 0) =0;
    
    protected:
        nthash::NtHash* hasher;

        ~Digester(){
            delete hasher;
        }
};

}

#endif
