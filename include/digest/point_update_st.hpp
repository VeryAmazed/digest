#ifndef SEG_TREE_HPP
#define SEG_TREE_HPP

#include <algorithm>
#include <vector>
#include <utility>
#include <stdint.h>


namespace segtree{
// Based on a template taken from USACO.guide and then modified by me (for competitive programming), and now modified again (for this)
// https://usaco.guide/gold/PURS?lang=cpp
// https://codeforces.com/blog/entry/18051 (USACO.guide was probably heavily inspired by this)
/** A data structure that can answer point update & range minimum queries. */

struct SegTree {
    // Max possible value for uint64_t
    const uint64_t uint64_t_MAX = 0xFFFFFFFFFFFFFFFF;
    
    // Default value, loses to everything
    const std::pair<uint64_t, size_t> DEFAULT = std::make_pair(uint64_t_MAX, 0);

    // array representation of complete binary tree
	std::vector<std::pair<uint64_t, size_t>> segtree;
	
    // number of leaves
    int len;

    SegTree(int len) : len(len){
        segtree = std::vector<std::pair<uint64_t, size_t>>(len * 2, DEFAULT);
    }

    SegTree(const SegTree& copy) : len(copy.len) {
        this->segtree = std::vector<std::pair<uint64_t, size_t>>(copy.segtree);
    }

    SegTree& operator=(const SegTree& copy){
        this->len = copy.len;
        segtree.assign(copy.segtree.begin(), copy.segtree.end());
        return *this;
    }

    ~SegTree(){
    }

    /** 
     * 
     * @param a the left value to be considered
     * @param b the right value to be considered
     * 
     * @return the object a and b with the minimum hash value (uint64_t), ties broken with the larger index
     *         ties are broken this way because we need default to always lose, but size_t doesn't have a well
     *         defined maximum value, so this allows us to make the default size_t value 0
	 */
	std::pair<uint64_t, size_t> comb(std::pair<uint64_t, size_t> a, std::pair<uint64_t, size_t> b) { 
		if(a.first < b.first){
            return a;
        }else if(a.first > b.first){
            return b;
        }else{
            if(a.second >= b.second){
                return a;
            }else{
                return b;
            }
        }
	}

	// look at ASSERTS to see how you should be indexing things
	/** 
     * @brief sets the value at ind to val. 
     * 
     * @param ind the index within the original sequence to be considered
     * @param val the value to set to
     */
	void set(size_t ind, std::pair<uint64_t, size_t> val) {
		// assert(0 <= ind && ind < len);
		ind += len;
		segtree[ind] = val;
		for (; ind > 1; ind /= 2) {
			segtree[ind >> 1] = comb(segtree[ind], segtree[ind ^ 1]);
		}
	}

	/** 
     * queries the range [start, end) 
     * @param start the 0-indexed indice indicating the inclusive left point of the range
     * @param end the 0-indexed indice indicating the exclusive right point of the range
     */
	// std::pair<uint64_t, size_t> query(size_t start, size_t end) {
	// 	// assert(0 <= start && start < len && 0 < end && end <= len);
	// 	std::pair<uint64_t, size_t> minAm = DEFAULT;
	// 	for (start += len, end += len; start < end; start /= 2, end /= 2) {
	// 		if ((start & 1) != 0) { minAm = comb(minAm, (*segtree)[start++]); }
	// 		if ((end & 1) != 0) { minAm = comb(minAm, (*segtree)[--end]); }
	// 	}
	// 	return minAm;
	// }
};

}
#endif
