#include "digest/thread_out.hpp"

namespace thread_out
{

void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint64_t mod, uint64_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_kmers = (int)len - (int)start - (int)k + 1;
        if(k < 4 ||start >= len || num_kmers < 0 || num_kmers < (int)thread_count){
            throw BadThreadOutParams();
        }
        unsigned kmers_per_thread = num_kmers/thread_count;
        unsigned extras = num_kmers % thread_count;
        vec.assign(thread_count, std::vector<size_t>());
        std::vector<std::thread> thread_vector;
        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            unsigned assigned_kmer_am = kmers_per_thread;
            if(extras > 0){
                ++assigned_kmer_am;
                extras--;
            }
            digest::ModMin dig(seq, assigned_kmer_am+k-1, k, mod, congruence, ind, minimized_h);
            thread_vector.push_back(std::thread(worker_roll, std::ref(dig), std::ref(vec[i]), std::ref(assigned_kmer_am)));
            ind += assigned_kmer_am;
        }
        for(auto& t: thread_vector)
        {
            t.join();
        }
    }
void thread_mod(unsigned thread_count, std::vector<std::vector<size_t>>& vec, 
    const std::string& seq, unsigned k, uint64_t mod, uint64_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_mod(thread_count, vec, seq.c_str(), seq.size(), k, mod, congruence, start, minimized_h);
    }


// small issue, thread args have to be copied in
void worker_roll(digest::Digester& dig, std::vector<size_t>& vec, unsigned& amount){
    dig.roll_minimizer(amount, vec);
}

}



