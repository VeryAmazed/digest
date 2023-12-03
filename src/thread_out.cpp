#include "digest/thread_out.hpp"

namespace thread_out
{

void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t mod, uint32_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_kmers = (int)len - (int)start - (int)k + 1;
        if(k < 4 ||start >= len || num_kmers < 0 || (unsigned)num_kmers < thread_count){
            throw BadThreadOutParams();
        }
        unsigned kmers_per_thread = num_kmers/thread_count;
        unsigned extras = num_kmers % thread_count;
        vec.assign(thread_count, std::vector<size_t>());
        std::vector<std::thread> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_kmer_am = kmers_per_thread;
            if(extras > 0){
                ++(assigned_kmer_am);
                extras--;
            }

            thread_vector.emplace_back(std::thread(thread_mod_roll, std::ref(vec[i]), 
                seq, ind, k, mod, congruence, minimized_h, assigned_kmer_am));

            ind += assigned_kmer_am;
        }
        for(auto& t: thread_vector)
        {
            t.join();
        }
    }

void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t mod, uint32_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_mod(thread_count, vec, seq.c_str(), seq.size(), k, mod, congruence, start, minimized_h);
    }

void thread_mod_roll(std::vector<size_t>& vec, const char* seq, 
    size_t ind, unsigned k, uint32_t mod, uint32_t congruence, 
    digest::MinimizedHashType minimized_h, unsigned assigned_kmer_am){
        digest::ModMin dig(seq, ind + assigned_kmer_am + k -1, k, mod, congruence, ind, minimized_h);
        dig.roll_minimizer(assigned_kmer_am, vec);
    }

}



