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
        UM_Digester(const std::string& seq, unsigned k, uint64_t mod, uint64_t congruence, size_t pos = 0, unsigned minimized_h = 0)
        :  Digester(seq, k, pos, minimized_h), mod(mod), congruence(congruence)
        {
            if(congruence >= mod){
                throw BadModException();
            }
        }

        UM_Digester(const char* seq, size_t len, unsigned k, uint64_t mod, uint64_t congruence, size_t pos = 0, unsigned minimized_h = 0)
        :  Digester(seq, len, k, pos, minimized_h), mod(mod), congruence(congruence)
        {
            if(congruence >= mod){
                throw BadModException();
            }
        }

    bool roll_next_minimizer() override;     
    
    private:
        uint64_t mod;
        uint64_t congruence;

};

}

#endif