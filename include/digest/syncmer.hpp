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
        Syncmer(const char* seq, size_t len, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON)
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
        Syncmer(const std::string& seq, unsigned k, unsigned large_wind_kmer_am, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON) :
            Syncmer(seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h)
        {}

        Syncmer(const Syncmer& copy) : WindowMin(copy){}

        Syncmer& operator=(const Syncmer& copy){
            WindowMin::operator=(copy);
            return *this;
        }

        /**
         * @brief adds up to amount of positions of syncmers into vec, here a large window is considered a syncmer if the smallest hash in the large window is at the leftmost or rightmost position 
         *        Time Complexity: O(log(large_wind_kmer_am)) per k-mer tested
         * 
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<size_t>& vec) override;

        /**
         * @brief adds up to amount of positions of syncmers into vec, here a large window is considered a syncmer if the smallest hash in the large window is at the leftmost or rightmost position 
         *        Time Complexity: O(log(large_wind_kmer_am)) per k-mer tested
         *        Because the size of syncmers may be irregular if there are k-mers with invalid characters within the large window that we skip over, this function adds pairs with the
         *        first element being the position of the first k-mer in the large window, and the second element being the position of the last k-mer in the large window, so if the distance
         *        between these two elements is not large_wind_kmer_am-1, then you know the window size is irregular, and you may need to do some processing on your end to get a syncmer of standard
         *        size that doesn't contain any k-mers with non-ACTG characters in it
         * 
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<std::pair<size_t, size_t>>& vec); 

    private:
        /**
         * @brief helper function that checks if a large_window is a syncmer for the function that takes a reference to a vector of size_t as its parameter
         * 
         * @param vec 
         */
        void check1(std::vector<size_t>& vec);

        /**
         * @brief helper function that checks if a large_window is a syncmer for the function that takes a reference to a vector of pairs of size_t's as its parameter
         * 
         * @param vec 
         */
        void check2(std::vector<std::pair<size_t, size_t>>& vec);
};

}

#endif