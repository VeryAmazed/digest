#ifndef SYNC_HPP
#define SYNC_HPP
#include "digester.hpp"
#include "window_minimizer.hpp"

namespace digest{

// number of k-mers to be considered in the large window
template <class P, class T>
class Syncmer : public WindowMin<P, T> {
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
         * @brief 
         * 
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<uint32_t>& vec) override;

        /**
         * @brief adds up to amount of positions and hashes of syncmers into vec, here a large window is considered a syncmer if the smallest hash in the large window is at the leftmost or rightmost position 
         *       
         * @param amount 
         * @param vec 
         */
        void roll_minimizer(unsigned amount, std::vector<std::pair<uint32_t, uint32_t>>& vec); 

    private:
        /**
         * @brief helper function which handles adding the next hash into the data structure
         * 
         */
		    void roll_ds_sync(std::vector<uint32_t>& vec);
        
        /**
         * @brief helper function which handles adding the next hash into the data structure
         * 
         */
        void roll_ds_sync(std::vector<std::pair<uint32_t, uint32_t>>& vec);
};

}

#include "syncmer.tpp"
#endif
