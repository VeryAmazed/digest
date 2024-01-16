#ifndef SYNC_HPP
#define SYNC_HPP
#include "digester.hpp"
#include "window_minimizer.hpp"

namespace digest{

// number of k-mers to be considered in the large window
template <class T>
class Syncmer : public WindowMin<T> {
    public:
      /**
         * 
         * @param seq 
         * @param len
         * @param k 
         * @param large_window 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        Syncmer(const char* seq, size_t len, unsigned k, unsigned large_window, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON)
        :  WindowMin<T>(seq, len, k, large_window, start, minimized_h)
        {}

        /**
         * 
         * @param seq 
         * @param k 
         * @param large_window 
         * @param start
         * @param minimized_h 
         * 
         * @throws BadWindowException Thrown when congruence is greater or equal to mod
         */
        Syncmer(const std::string& seq, unsigned k, unsigned large_window, size_t start = 0, MinimizedHashType minimized_h = MinimizedHashType::CANON) :
            Syncmer(seq.c_str(), seq.size(), k, large_window, start, minimized_h)
        {}

        /**
         * @brief adds up to amount of positions of syncmers into vec, here a large window is considered a syncmer if the smallest hash in the large window is at the leftmost or rightmost position 
         * 
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<uint32_t>& vec) override;

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
        void roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec); 

    private:
		void fill_st(std::vector<uint32_t>& vec);
		void fill_st(std::vector<std::pair<uint32_t, uint32_t>>& vec);
};

}

#include "syncmer.tpp"
#endif
