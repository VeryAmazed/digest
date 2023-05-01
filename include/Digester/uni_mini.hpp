#ifndef UNI_MINI_HPP
#define UNI_MINI_HPP

#include "digester.hpp"

namespace digest{

class BadModException : public std::exception
{
	const char * what () const throw ()
    {
    	return "mod must be greater than congruence.";
    }
};

class ModMin : public Digester{
    public:
        /**
         * 
         * @param seq 
         * @param len
         * @param k 
         * @param mod Mod space to be used to calculate universal minimizers
         * @param congruence value we want minimizer hashes to be congruent to in the mod space
         * @param start
         * @param minimized_h 
         * 
         * @throws BadModException Thrown when congruence is greater or equal to mod
         */
        ModMin(const char* seq, size_t len, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t start = 0, unsigned minimized_h = 0)
        :  Digester(seq, len, k, start, minimized_h), mod(mod), congruence(congruence)
        {
            if(congruence >= mod){
                throw BadModException();
            }
        }

        /**
         * 
         * @param seq 
         * @param k 
         * @param mod Mod space to be used to calculate universal minimizers
         * @param congruence value we want minimizer hashes to be congruent to in the mod space
         * @param start
         * @param minimized_h 
         * 
         * @throws BadModException Thrown when congruence is greater or equal to mod
         */
        ModMin(const std::string& seq, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t start = 0, unsigned minimized_h = 0) :
            ModMin(seq.c_str(), seq.size(), k, mod, congruence, start, minimized_h)
        {}
        

        /**
         * Copy Constructor
         * 
         * @param copy, ModMin object you want to copy from 
         */
        ModMin(const ModMin& copy) : Digester(copy), mod(copy.mod), congruence(copy.congruence)
        {}

        ModMin& operator=(const ModMin& copy){
            this->mod = copy.mod;
            this->congruence = copy.congruence;
            Digester::operator=(copy);
            return *this;
        }

        std::vector<size_t> roll_minimizer(unsigned amount) override;    

        /**
         * 
         * @return uint64_t, the mod space being used
         */
        uint64_t get_mod(){
            return mod;
        }

        /**
         * 
         * @return uint64_t, the value the minimized hash must be congruent to
         */
        uint64_t get_congruence(){
            return congruence;
        } 
    
    private:
        uint64_t mod;
        uint64_t congruence;

};

}

#endif