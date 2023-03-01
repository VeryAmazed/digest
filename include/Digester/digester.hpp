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
        /*
        Only supports std::string for now
        */
        Digester(const std::string& seq, unsigned k, size_t pos = 0, bool canonicalized){
            hasher = new nthash::NtHash(seq, 1, k, pos);
            this->canonicalized = canonicalized;
        }

        virtual bool roll() = 0;

        virtual uint64_t get_minimizer() = 0;

        virtual void change_seq(const std::string& new_seq, size_t new_pos = 0) =0;

        virtual ~Digester(){
            delete hasher;
        }
        
    protected:
        nthash::NtHash* hasher;

    private:
        bool canonicalized;
};

}

#endif
