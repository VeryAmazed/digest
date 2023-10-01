#pragma once

#include "nthash/internal.hpp"
#include "nthash/nthash.hpp"

namespace {
/**
 * Check the current k-mer for non ACGTU's
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return `true` if any of the first k characters is not an ACGTU, `false`
 * otherwise
 */
inline bool
is_invalid_kmer(const char* seq, unsigned k, size_t& pos_n);

/**
 * Generate the forward-strand hash value of the first k-mer in the sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of k-mer_0
 */
inline uint64_t
base_forward_hash(const char* seq, unsigned k);

/**
 * Perform a roll operation on the forward strand by removing char_out and
 * including char_in.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled forward hash value
 */
inline uint64_t
next_forward_hash(uint64_t fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in);

/**
 * Perform a roll back operation on the forward strand.
 * @param fh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Forward hash value rolled back
 */
inline uint64_t
prev_forward_hash(uint64_t fh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in);

/**
 * Generate a hash value for the reverse-complement of the first k-mer in the
 * sequence.
 * @param seq C array containing the sequence's characters
 * @param k k-mer size
 * @return Hash value of the reverse-complement of k-mer_0
 */
inline uint64_t
base_reverse_hash(const char* seq, unsigned k);

/**
 * Perform a roll operation on the reverse-complement by removing char_out and
 * including char_in.
 * @param rh_val Previous reverse-complement hash value computed for the
 * sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Rolled hash value for the reverse-complement
 */
inline uint64_t
next_reverse_hash(uint64_t rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in);

/**
 * Perform a roll back operation on the reverse strand.
 * @param rh_val Previous hash value computed for the sequence
 * @param k k-mer size
 * @param char_out Character to be removed
 * @param char_in Character to be included
 * @return Reverse hash value rolled back
 */
inline uint64_t
prev_reverse_hash(uint64_t rh_val,
                  unsigned k,
                  unsigned char char_out,
                  unsigned char char_in);
} // namespace
