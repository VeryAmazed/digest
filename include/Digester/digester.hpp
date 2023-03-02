#ifndef DIGESTER_HPP
#define DIGESTER_HPP

#include "../ntHash/nthash.hpp"
#include <stdexcept>
#include <deque>

namespace digest{

// Only supports characters in DNA
class UnknownCharException : public std::exception
{
	const char * what () const throw ()
    {
    	return "Sequence contains character that is not A, C, T, or G (upper or lower case)";
    }
};

class NotRolledException : public std::exception
{
	const char * what () const throw ()
    {
    	return "Must call roll_one() or roll_next_minimizer() first.";
    }
};

// Only supports characters in DNA, upper or lower case
class Digester{
    public:
        
        Digester(const std::string& seq, unsigned k, size_t start = 0, uint8_t minimized_h = 0) 
            : seq(seq.data()), len(seq.size()), pos(start), start(start), end(start+k), k(k), minimized_h(minimized_h)
        {
           
        }

        Digester(const char* seq, size_t len, unsigned k, size_t start = 0, uint8_t minimized_h = 0) 
            : seq(seq), len(len), pos(start), start(start), end(start+k), k(k), minimized_h(minimized_h)
        {
            
        }

        void roll_one();

        virtual void roll_next_minimizer() = 0;

        size_t get_pos(){
            return pos;
        }

        uint64_t get_chash(){
            if(!rolled){
                throw NotRolledException();
            }
            return chash;
        }

        uint64_t get_fhash(){
            if(!rolled){
                throw NotRolledException();
            }
            return fhash;
        }

        uint64_t get_rhash(){
            if(!rolled){
                throw NotRolledException();
            }
            return rhash;
        }
        
        void new_seq(const std::string& seq, size_t pos){
            this->seq = seq.data();
            this->len = len;
            this->pos = pos;
            this->start = pos;
            this->end = pos+this->k;
            rolled = false;
        }

        void new_seq(const char* seq, size_t len, size_t pos){
            this->seq = seq;
            this->len = len;
            this->pos = pos;
            this->start = pos;
            this->end = pos+this->k;
            rolled = false;
        }

        void append_seq(const std::string& seq);

        void append_seq(const char* seq, size_t len);

        uint8_t get_minimized_h(){
            return minimized_h;
        }

        // to be added
        std::string get_string();

    protected:
        // sequence to be digested
        const char* seq;
        // length of seq
        size_t len;
        // pos within entirety of the sequence you are digesting, sequences that are appended by append_seq are counted as one sequence
        size_t pos;
        // internal index of the next character to be thrown out, junk if c_outs is not empty
        size_t start;
        // internal index of next character to be added
        size_t end;
        // canonical hash
        uint64_t chash;
        // forward hash
        uint64_t fhash;
        // reverse hash
        uint64_t rhash;
        // length of kmer
        unsigned k;
        // deque of characters to be thrown out from left to right
        std::deque<char> c_outs;

    private:
        /*
            Hash value to be minimized
            0 for canonical, 1 for forward, 2 for reverse
        */
        uint8_t minimized_h;
        // internal bool to track if rolled was called at least once
        bool rolled = false;
};

}

#endif
