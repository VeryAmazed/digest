#ifndef DIGESTER_HPP
#define DIGESTER_HPP

#include <nthash/nthash.hpp>
#include <nthash/kmer.hpp>
#include <stdexcept>
#include <deque>
#include <cctype>

namespace digest{

class BadConstructionException : public std::exception
{
	const char * what () const throw ()
    {
    	return "minimized_h must be either 0, 1, or 2, k must be greater than 3, start must be less than len";
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
         * @param seq char pointer poitning to the c-string of DNA sequence to be hashed.
         * @param len length of seq.
         * @param k k-mer size.
         * @param start 0-indexed position in seq to start hashing from. 
         * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
         * 
         * @throws BadConstructionException Thrown if k equals 0 or is greater than the length of the sequence, if minimized_h is not 0, 1, or 2,
         *      or if the starting position is not at least k-1 from the end of the string
         */
        Digester(const char* seq, size_t len, unsigned k, size_t start = 0, unsigned minimized_h = 0) 
            : seq(seq), len(len), offset(0), start(start), end(start+k), k(k), minimized_h(minimized_h), fhash(0), chash(0), rhash(0) {
                if(k < 4 ||start >= len || minimized_h > 2){
                    throw BadConstructionException();
                }
                this->c_outs = new std::deque<char>;
                init_hash();
            }
        
        /**
         * @param seq reference to std string of DNA sequence to be hashed.
         * @param k k-mer size. 
         * @param start 0-indexed position in seq to start hashing from. 
         * @param minimized_h hash to be minimized, 0 for canoncial, 1 for forward, 2 for reverse
         * 
         * @throws BadConstructionException Thrown if k equals 0 or is greater than the length of the sequence, if minimized_h is not 0, 1, or 2,
         *      or if the starting position is not at least k-1 from the end of the string
         */
        Digester(const std::string& seq, unsigned k, size_t start = 0, unsigned minimized_h = 0) :
            Digester(seq.c_str(), seq.size(), k, start, minimized_h) {}

        Digester(const Digester& copy){
            copyOver(copy);
            this->c_outs = new std::deque<char>(*(copy.c_outs));
        }

        Digester& operator=(const Digester& copy){
            copyOver(copy);
            (this->c_outs)->assign((copy.c_outs)->begin(), (copy.c_outs)->end());
            return *this;
        }

        virtual ~Digester(){
            delete c_outs;
        }
        
        /**
         * @return bool, true if values of the 3 hashes are meaningful, false otherwise, i.e. the object wasn't able to initialize with a valid hash or roll_one() was called when already at end of sequence
         */
        bool get_is_valid_hash(){
            return is_valid_hash;
        }

        unsigned get_k(){
            return k;
        }

        size_t get_len(){
            return len;
        }

        /**
         * @brief moves the internal pointer to the next valid k-mer, skipping over any k-mers that have contain a non ACTG character, and returns hash for that k-mer 
         *        Time Complexity: O(1)
         * 
         * @return bool, true if we were able generate a valid hash, false otherwise
         */
        bool roll_one();

        /**
         * @brief returns the positions, as defined by get_pos(), of minimizers up to the amount specified 
         * 
         * @param amount number of minimizers you want to generate
         * @param vec a reference to a vector of size_t's, the positions returned will go there
         */
        virtual void roll_minimizer(unsigned amount, std::vector<size_t>& vec) = 0;

        /**
         * @return current index of the first character of the current kmer that has been hashed
         *         strings that have been appended onto each other count as 1 big string, 
         *         i.e. if you first had a string of length 10 and then appended another string of length 20, and the index of the first character of the current
         *         k-mer is at index 4, 0-indexed, in the second string, then it will return 14
         */
        size_t get_pos(){
            return offset + start - c_outs->size();
        }

        uint64_t get_chash(){
            return chash;
        }

        uint64_t get_fhash(){
            return fhash;
        }

        uint64_t get_rhash(){
            return rhash;
        }

        /**
         * @brief replaces the current sequence with the new one, it's like starting over with a completely new string
         * 
         * @param seq char pointer to new sequence to be hashed
         * @param len length of the new sequence
         * @param start position in new sequence to start from
         * 
         * @throws BadConstructionException thrown if the starting position is greater than the length of the string
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
         * @brief replaces the current sequence with the new one, it's like starting over with a completely new string
         * 
         * @param seq std string reference to the new sequence to be hashed
         * @param start position in new sequence to start from
         * 
         * @throws BadConstructionException thrown if the starting position is greater than the length of the string
         */
        void new_seq(const std::string& seq, size_t pos){
            new_seq(seq.c_str(), seq.size(), pos);
        }

        /**
         * @brief simulates the appending of a new sequence to the end of the old sequence
         * The old string will no longer be stored, but the rolling hashes will be able to preceed as if the strings were appended
         * Can only be called when you've reached the end of the current string
         * i.e. if you're current sequence is ACTGAC, and you have reached the end of this sequence, and you call append_seq with the
         * sequence CCGGCCGG, then the minimizers you will get after calling append_seq plus the minimizers you got from
         * going through ACTGAC, will be equivalent to the minimizers you would have gotten from rolling across ACTGACCCGGCCGG
         * 
         * @param seq C string of DNA sequence to be appended
         * @param len length of the sequence
         * 
         * @throws NotRolledTillEndException Thrown when the internal iterator is not at the end of the current sequence
         */
        void append_seq(const char* seq, size_t len);

        /**
         * @brief simulates the appending of a new sequence to the end of the old sequence
         * The old string will no longer be stored, but the rolling hashes will be able to preceed as if the strings were appended
         * Can only be called when you've reached the end of the current string
         * i.e. if you're current sequence is ACTGAC, and you have reached the end of this sequence, and you call append_seq with the
         * sequence CCGGCCGG, then the minimizers you will get after calling append_seq plus the minimizers you got from
         * going through ACTGAC, will be equivalent to the minimizers you would have gotten from rolling across ACTGACCCGGCCGG
         * 
         * @param seq std string of DNA sequence to be appended
         * 
         * @throws NotRolledTillEndException Thrown when the internal iterator is not at the end of the current sequence
         */
        void append_seq(const std::string& seq);

        /**
         * @return unsigned, a number representing the hash you are minimizing, 0 for canoncial, 1 for forward, 2 for reverse 
         */
        unsigned get_minimized_h(){
            return minimized_h;
        }

        /**
         * @return const char* representation of the sequence
         */
        const char* get_sequence(){
            return seq;
        }
        
    protected:
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
         * Helper function
         * 
         * @param in char to be checked 
         * @return bool, true if in is an upper or lowercase ACTG character, false otherwise
         */
        bool is_ACTG(char in){
            in = toupper(in);
            if(in == 'A' || in == 'C' || in == 'T' || in == 'G'){
                return true;
            }
            return false;
        }

        /**
         * @brief Helper function that initializes the hash values at the first valid k-mer at or after start
         *        Sets is_valid_hash to be equal to its return value
         * 
         * @return bool, true on success, a valid hash is initialized, false otherwise
         */
        bool init_hash();


        // sequence to be digested, memory is owned by the user
        const char* seq;
        
        // length of seq
        size_t len;
        
        // the combined length of all the previous strings that have been appended together, not counting the current string
        size_t offset;
        
        // internal index of the next character to be thrown out, junk if c_outs is not empty
        size_t start;
        
        // internal index of next character to be added
        size_t end;
        
        // canonical hash of current k-mer
        uint64_t chash;
        
        // forward hash of current k-mer
        uint64_t fhash;
        
        // reverse hash of current k-mer
        uint64_t rhash;
        
        // length of kmer
        unsigned k;
        
        // deque of characters to be rolled out in the rolling hash from left to right, memory is owned by the object
        std::deque<char>* c_outs;
        
        //Hash value to be minimized, 0 for canonical, 1 for forward, 2 for reverse
        unsigned minimized_h;

        // bool representing whether the current hash is meaningful, i.e. corresponds to the k-mer at get_pos()
        bool is_valid_hash = false;

};

}

#endif
