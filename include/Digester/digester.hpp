#ifndef DIGESTER_HPP
#define DIGESTER_HPP

#include "../ntHash/nthash.hpp"

namespace digest{

class Digester{
    public:
        /*
        Only supports std::string for now
        */
        Digester(const std::string& seq, unsigned k, size_t pos = 0, bool canonicalized = true){
            hasher = new nthash::NtHash(seq, 1, k, pos);
            this->canonicalized = canonicalized;
            this->seq = seq;
        }

        // take you to next minimizer, first item in pair denotes success/failure, second in position in seq
        virtual std::pair<bool, size_t> roll() = 0;

        // only after you call roll, gives you minimizer string
        virtual std::string get_min_string() = 0;

        // might just be for testing
        // only after you call roll, gives you minimizer hash
        virtual uint64_t get_min_hash() =0;

        virtual void change_seq(const std::string& new_seq, size_t new_pos) =0;

        virtual ~Digester(){
            delete hasher;
        }
        
    protected:
        nthash::NtHash* hasher;
        std::string seq;

    private:
        bool canonicalized;
};

}

#endif
