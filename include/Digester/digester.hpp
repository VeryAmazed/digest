#ifndef DIGESTER_HPP
#define DIGESTER_HPP

#include "../ntHash/nthash.hpp"
#include <stdexcept>
#include <deque>
#include <cctype>

namespace digest{

// Only supports characters in DNA and N

class BadConstructionException : public std::exception
{
	const char * what () const throw ()
    {
    	return "minimized_h must be either 0, 1, or 2, k cannot be 0, start must be less than len";
    }
};

class NotRolledTillEndException : public std::exception
{
	const char * what () const throw ()
    {
    	return "Iterator must be at the end of the current sequence before appending a new one.";
    }
};

// Only supports characters in DNA and N, upper or lower case
class Digester{
    public:

        /**
         * Constructor.
         * @param seq C string of DNA sequence to be hashed.
         * @param len Length of seq.
         * @param k K-mer size.
         * @param start 0-indexed position in seq to start hashing from. 
         * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
         * 
         * @throws BadConstructionException Thrown if k equals 0 or is greater than the length of the sequence, if minimized_h is not 0, 1, or 2,
         *      or if the starting position is not at least k-1 from the end of the string
         */
        Digester(const char* seq, size_t len, unsigned k, size_t start = 0, unsigned minimized_h = 0) 
            : seq(seq), len(len), offset(0), start(start), end(start+k), k(k), minimized_h(minimized_h) {
                fhash = 0;
                chash = 0;
                rhash = 0;
                if(k == 0 ||start >= len || minimized_h > 2){
                    throw BadConstructionException();
                }
                this->c_outs = new std::deque<char>;
                init_hash();
            }
        
        /**
         * Constructor.
         * @param seq std string of DNA sequence to be hashed.
         * @param k K-mer size. 
         * @param start 0-indexed position in seq to start hashing from. 
         * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
         * 
         * @throws BadConstructionException Thrown if k equals 0 or is greater than the length of the sequence, if minimized_h is not 0, 1, or 2,
         *      or if the starting position is not at least k-1 from the end of the string
         */
        Digester(const std::string& seq, unsigned k, size_t start = 0, unsigned minimized_h = 0) :
            Digester(seq.c_str(), seq.size(), k, start, minimized_h) {}

        /**
         * Helper function
         * 
         * @param copy, Digester object you want to copy from 
         */
        void copyOver(const Digester& copy){
            this->seq = copy.seq;
            this->len = copy.len;
            this->k = copy.k;
            this->offset = copy.offset;
            this->start = copy.start;
            this->end = copy.end;
            this->is_valid_hash = copy.is_valid_hash;
            this->minimized_h = copy.minimized_h;
            if(this->is_valid_hash){
                this->chash = copy.chash;
                this->rhash = copy.rhash;
                this->fhash = copy.fhash;
            }
        }
        /**
         * Copy Constructor
         * 
         * @param copy, Digester object you want to copy from 
         */
        Digester(const Digester& copy){
            copyOver(copy);
            this->c_outs = new std::deque<char>(*(copy.c_outs));
        }
        /**
         * assignment operator override
         * 
         * @param copy, Digester object you want to copy from 
         */
        Digester& operator=(const Digester& copy){
            copyOver(copy);
            (this->c_outs)->assign((copy.c_outs)->begin(), (copy.c_outs)->end());
            return *this;
        }

        virtual ~Digester(){
            delete c_outs;
        }
        
        /**
         * @return bool, true if roll_one(), roll_next_minimizer() or roll_next_n_minis() has been called at least once, false otherwise
         * 
         */
        bool get_is_valid_hash(){
            return is_valid_hash;
        }
        
        /**
         * 
         * @return unsigned value of k
         */
        unsigned get_k(){
            return k;
        }

        /**
         * 
         * @return size_t, length of the sequence
         */
        size_t get_len(){
            return len;
        }

        /**
         * roll the hash 1offsetition to the right or construcuts the initial hash on first call 
         * 
         * @throws std::out_of_range if the end of the string has already been reached
         */
        bool roll_one();

        /**
         * 
         * @param amount number of minimizers you want to generate
         * @return std::vector<size_t> vector filled with the positions of the minimizers up to the amount
         */
        virtual std::vector<size_t> roll_minimizer(unsigned amount) = 0;

        /**
         * 
         * @return current index of the first character of the kmer that has been hashed
         *         strings that have been appended onto each other count as 1 big string
         */
        size_t get_pos(){
            return offset + start - c_outs->size();
        }

        /**
         * 
         * @return the canonicalized hash of the current k-mer
         */
        uint64_t get_chash(){
            return chash;
        }

        /**
         * 
         * @return the forward hash of the current k-mer
         */
        uint64_t get_fhash(){
            return fhash;
        }

        /**
         * 
         * @return the reverse hash of the current k-mer
         */
        uint64_t get_rhash(){
            return rhash;
        }

        /**
         * 
         * @param seq new sequence to be hashed
         * @param len length of the new sequence
         * @paramoffset newoffsetition to start from
         * 
         * @throws BadConstructionException Thrown if k is greater than the length of the sequence,
         *      or if the startingoffsetition is not at least k-1 from the end of the string
         */
        void new_seq(const char* seq, size_t len, size_t start){
            this->seq = seq;
            this->len = len;
            this->offset = 0;
            this->start = start;
            this->end = start+this->k;
            is_valid_hash = false;
            if(start >= len){
                throw BadConstructionException();
            }
            init_hash();
        }

        /**
         * 
         * @param seq new sequence to be hashed
         * @paramoffset newoffsetition to start from
         * 
         * @throws BadConstructionException Thrown if k is greater than the length of the sequence,
         *      or if the startingoffsetition is not at least k-1 from the end of the string
         */
        void new_seq(const std::string& seq, size_t pos){
            new_seq(seq.c_str(), seq.size(), pos);
        }

        /**
         * Simulates the appending of a new sequence to the end of the old sequence
         * The old string will no longer be stored, but the rolling hashes will be able to preceed as if the strings were appended
         * 
         * @param seq C string of DNA sequence to be appended
         * @param len length of the sequence
         * 
         * @throws NotRolledTillEndException Thrown when the internal iterator is not at the end of the current sequence
         */
        void append_seq(const char* seq, size_t len);

        /**
         * Simulates the appending of a new sequence to the end of the old sequence
         * The old string will no longer be stored, but the rolling hashes will be able to preceed as if the strings were appended
         * 
         * @param seq std string of DNA sequence to be appended
         * 
         * @throws NotRolledTillEndException Thrown when the internal iterator is not at the end of the current sequence
         */
        void append_seq(const std::string& seq);

        unsigned get_minimized_h(){
            return minimized_h;
        }

        bool is_ACTG(char in){
            in = toupper(in);
            if(in == 'A' || in == 'C' || in == 'T' || in == 'G'){
                return true;
            }
            return false;
        }

        // Helper function that initializes the hash values
        bool init_hash();

        /**
         * 
         * @return const char* representation of the sequence
         */
        const char* get_sequence(){
            return seq;
        }

        
        /**
         * 
         * @return std::string of the current k-mer
         */
        // std::string get_string();
        
    protected:
        // sequence to be digested
        const char* seq;
        
        // length of seq
        size_t len;
        
        // the combined length of all the previous strings that have been appended together, not counting the current string
        size_t offset;
        
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
        std::deque<char>* c_outs;
        
        /*
            Hash value to be minimized
            0 for canonical, 1 for forward, 2 for reverse
        */
        unsigned minimized_h;

        // internal bool to track if rolled was called at least once
        bool is_valid_hash = false;
};

}

#endif
