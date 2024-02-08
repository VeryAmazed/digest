#ifndef MOD_MINI_HPP
#define MOD_MINI_HPP

#include "digester.hpp"
#include <cassert>
#include <cstdint>

namespace digest{

class BadModException : public std::exception
{
	const char * what () const throw ()
    {
    	return "mod must be greater than congruence.";
    }
};
template <class P>
class ModMin : public Digester<P>{
    public:
        /**
         * 
         * @param seq 
         * @param len
         * @param k 
         * @param mod mod space to be used to calculate universal minimizers
         * @param congruence value we want minimizer hashes to be congruent to in the mod space
         * @param start
         * @param minimized_h 
         * 
         * @throws BadModException Thrown when congruence is greater or equal to mod
         */
        ModMin(const char* seq, size_t len, unsigned k, uint32_t mod, uint32_t congruence = 0, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON)
        :  Digester(seq, len, k, start, minimized_h), mod(mod), congruence(congruence)
        {
            if (congruence >= mod){
                throw BadModException();
            }
        }

        /**
         * 
         * @param seq 
         * @param k 
         * @param mod mod space to be used to calculate universal minimizers
         * @param congruence value we want minimizer hashes to be congruent to in the mod space
         * @param start
         * @param minimized_h 
         * 
         * @throws BadModException Thrown when congruence is greater or equal to mod
         */
        ModMin(const std::string& seq, unsigned k, uint32_t mod, uint32_t congruence = 0, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON) :
            ModMin(seq.c_str(), seq.size(), k, mod, congruence, start, minimized_h)
        {}
        
        /**
         * @brief adds up to amount of positions of minimizers into vec, here a k-mer is considered a minimizer if its hash is congruent to congruence in the mod space 
         *        Time Complexity: O(1) per k-mer tested
         * 
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<uint32_t>& vec) override;    

         /**
         * @brief adds up to amount of positions and hashes of minimizers into vec, here a k-mer is considered a minimizer if its hash is congruent to congruence in the mod space 
         *        Time Complexity: O(1) per k-mer tested
         * 
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec) override;    


        /**
         * @return uint32_t, the mod space being used
         */
        uint32_t get_mod(){
            return mod;
        }

        /**
         * @return uint32_t, the value the minimized hash must be congruent to
         */
        uint32_t get_congruence(){
            return congruence;
        } 
    
    private:
        uint32_t mod;
        uint32_t congruence;

};

}

#endif
