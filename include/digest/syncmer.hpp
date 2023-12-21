#ifndef SYNC_HPP
#define SYNC_HPP
#include "digester.hpp"
#include "window_minimizer.hpp"

namespace digest{

// number of k-mers to be considered in the large window
template <uint32_t large_window>
class Syncmer : public WindowMin<large_window> {
    public:
      /**
         * 
         * @param seq 
         * @param len
         * @param k 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        Syncmer(const char* seq, size_t len, unsigned k, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON)
        :  WindowMin<large_window>(seq, len, k, start, minimized_h)
        {}

        /**
         * 
         * @param seq 
         * @param k 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        Syncmer(const std::string& seq, unsigned k, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON) :
            Syncmer(seq.c_str(), seq.size(), k, start, minimized_h)
        {}

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

#include "syncmer.tpp"
#endif
