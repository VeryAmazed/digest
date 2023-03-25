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

class UM_Digester : public Digester{
    public:

        /**
         * 
         * @param seq 
         * @param k 
         * @param mod Mod space to be used to calculate universal minimizers
         * @param congruence value we want minimizer hashes to be congruent to in the mod space
         * @param pos 
         * @param minimized_h 
         * 
         * @throws BadModException Thrown when congruence is greater or equal to mod
         */
        UM_Digester(const std::string& seq, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t pos = 0, unsigned minimized_h = 0)
        :  Digester(seq, k, pos, minimized_h), mod(mod), congruence(congruence)
        {
            if(congruence >= mod){
                throw BadModException();
            }
        }

        /**
         * 
         * @param seq 
         * @param len
         * @param k 
         * @param mod Mod space to be used to calculate universal minimizers
         * @param congruence value we want minimizer hashes to be congruent to in the mod space
         * @param pos 
         * @param minimized_h 
         * 
         * @throws BadModException Thrown when congruence is greater or equal to mod
         */
        UM_Digester(const char* seq, size_t len, unsigned k, uint64_t mod, uint64_t congruence = 0, size_t pos = 0, unsigned minimized_h = 0)
        :  Digester(seq, len, k, pos, minimized_h), mod(mod), congruence(congruence)
        {
            if(congruence >= mod){
                throw BadModException();
            }
        }

        bool roll_next_minimizer() override;    

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