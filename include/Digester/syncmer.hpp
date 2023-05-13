#ifndef SYNC_HPP
#define SYNC_HPP
#include "digester.hpp"
#include "window_minimizer.hpp"

namespace digest{

class Syncmer : public WindowMin{
    public:
      /**
         * 
         * @param seq 
         * @param len
         * @param k 
         * @param large_wind_kmer_am number of kmers to be considered at a time
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        Syncmer(const char* seq, size_t len, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, unsigned minimized_h = 0)
        :  WindowMin(seq, len, k, large_wind_kmer_am, start, minimized_h)
        {}

        /**
         * 
         * @param seq 
         * @param k 
         * @param large_wind_kmer_am number of kmers to be considered at a time
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        Syncmer(const std::string& seq, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, unsigned minimized_h = 0) :
            Syncmer(seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h)
        {}

        Syncmer(const Syncmer& copy) : WindowMin(copy){}

        Syncmer& operator=(const Syncmer& copy){
            WindowMin::operator=(copy);
            return *this;
        }

        void roll_minimizer(unsigned amount, std::vector<size_t>& vec) override;

        void roll_minimizer(unsigned amount, std::vector<std::pair<size_t, size_t>>& vec); 

    private:
        void check1(std::vector<size_t>& vec);

        void check2(std::vector<std::pair<size_t, size_t>>& vec);
};

}

#endif