#include "digest/thread_out.hpp"
#include <future>

namespace digest::thread_out
{
template<digest::BadCharPolicy P>
void thread_mod(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t mod, uint32_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_kmers = (int)len - (int)start - (int)k + 1;
        if(k < 4 ||start >= len || num_kmers < 0 || (unsigned)num_kmers < thread_count){
            throw BadThreadOutParams();
        }
        unsigned kmers_per_thread = num_kmers/thread_count;
        unsigned extras = num_kmers % thread_count;
		vec.reserve(thread_count);
        std::vector<std::future<std::vector<uint32_t>>> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_kmer_am = kmers_per_thread;
            if(extras > 0){
                ++(assigned_kmer_am);
                extras--;
            }

            thread_vector.emplace_back(std::async(thread_mod_roll1<P>,
                seq, ind, k, mod, congruence, minimized_h, assigned_kmer_am));

            ind += assigned_kmer_am;
        }
        for(auto& t: thread_vector)
        {
            vec.emplace_back(t.get());
        }
    }

template<digest::BadCharPolicy P>
void thread_mod(unsigned thread_count, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t mod, uint32_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_kmers = (int)len - (int)start - (int)k + 1;
        if(k < 4 ||start >= len || num_kmers < 0 || (unsigned)num_kmers < thread_count){
            throw BadThreadOutParams();
        }
        unsigned kmers_per_thread = num_kmers/thread_count;
        unsigned extras = num_kmers % thread_count;
        vec.reserve(thread_count);
        std::vector<std::future<std::vector<std::pair<uint32_t,uint32_t>>>> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_kmer_am = kmers_per_thread;
            if(extras > 0){
                ++(assigned_kmer_am);
                extras--;
            }

            thread_vector.emplace_back(std::async(thread_mod_roll2<P>,
                seq, ind, k, mod, congruence, minimized_h, assigned_kmer_am));

            ind += assigned_kmer_am;
        }
        for(auto& t: thread_vector)
        {
			vec.emplace_back(t.get());
        }
    }

template<digest::BadCharPolicy P>
void thread_mod(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t mod, uint32_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_mod<P>(thread_count, vec, seq.c_str(), seq.size(), k, mod, congruence, start, minimized_h);
    }

template<digest::BadCharPolicy P>
void thread_mod(unsigned thread_count, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& vec, 
    const std::string& seq, unsigned k, uint32_t mod, uint32_t congruence, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_mod<P>(thread_count, vec, seq.c_str(), seq.size(), k, mod, congruence, start, minimized_h);
    }

template<digest::BadCharPolicy P>
std::vector<uint32_t> thread_mod_roll1(const char* seq, 
    size_t ind, unsigned k, uint32_t mod, uint32_t congruence, 
    digest::MinimizedHashType minimized_h, unsigned assigned_kmer_am){
		std::vector<uint32_t> out;
        digest::ModMin<P> dig(seq, ind + assigned_kmer_am + k -1, k, mod, congruence, ind, minimized_h);
        dig.roll_minimizer(assigned_kmer_am, out);
		return out;
    }

template<digest::BadCharPolicy P>
std::vector<std::pair<uint32_t,uint32_t>> thread_mod_roll2(const char* seq, 
    size_t ind, unsigned k, uint32_t mod, uint32_t congruence, 
    digest::MinimizedHashType minimized_h, unsigned assigned_kmer_am){
		std::vector<std::pair<uint32_t,uint32_t>> out;
        digest::ModMin<P> dig(seq, ind + assigned_kmer_am + k -1, k, mod, congruence, ind, minimized_h);
        dig.roll_minimizer(assigned_kmer_am, out);
		return out;
    }

template <digest::BadCharPolicy P, class T>
void thread_wind(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_lwinds = (int)len - (int)start - (int)(k+large_wind_kmer_am)+2;
        if(large_wind_kmer_am == 0 || k < 4 ||start >= len || num_lwinds < 0 || (unsigned)num_lwinds < thread_count){
            throw BadThreadOutParams();
        }
        unsigned lwinds_per_thread = num_lwinds/thread_count;
        unsigned extras = num_lwinds % thread_count;
        vec.reserve(thread_count);
        std::vector<std::future<std::vector<uint32_t>>> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_lwind_am = lwinds_per_thread;
            if(extras > 0){
                ++(assigned_lwind_am);
                extras--;
            }

            thread_vector.emplace_back(std::async(thread_wind_roll1<P, T>,
                seq, ind, k, large_wind_kmer_am, minimized_h, assigned_lwind_am));

            ind += assigned_lwind_am;
        }
        for(auto& t: thread_vector)
        {
			vec.emplace_back(t.get());
        }
        
        // handle duplicates
        // the only possible place for a duplicate is for the last element
        // of vec[i] to equal the first value of vec[i+1] due to the fact
        // that thread_i+1 can't know the last minimizer of thread_i
        for(unsigned i = 0; i < thread_count-1; i++){
            int last = (int)vec[i].size() - 1;
            if(vec[i][last] == vec[i+1][0]){
                vec[i].pop_back();
            }
        }
    }

template <digest::BadCharPolicy P, class T>
void thread_wind(unsigned thread_count, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_lwinds = (int)len - (int)start - (int)(k+large_wind_kmer_am)+2;
        if(large_wind_kmer_am == 0 || k < 4 ||start >= len || num_lwinds < 0 || (unsigned)num_lwinds < thread_count){
            throw BadThreadOutParams();
        }
        unsigned lwinds_per_thread = num_lwinds/thread_count;
        unsigned extras = num_lwinds % thread_count;
        vec.reserve(thread_count);
        std::vector<std::future<std::pair<uint32_t, uint32_t>>> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_lwind_am = lwinds_per_thread;
            if(extras > 0){
                ++(assigned_lwind_am);
                extras--;
            }

            thread_vector.emplace_back(std::async(thread_wind_roll2<P, T>,
                seq, ind, k, large_wind_kmer_am, minimized_h, assigned_lwind_am));

            ind += assigned_lwind_am;
        }
        for(auto& t: thread_vector)
        {
			vec.emplace_back(t.get());
        }
        
        // handle duplicates
        // the only possible place for a duplicate is for the last element
        // of vec[i] to equal the first value of vec[i+1] due to the fact
        // that thread_i+1 can't know the last minimizer of thread_i
        for(unsigned i = 0; i < thread_count-1; i++){
            int last = (int)vec[i].size() - 1;
            if(vec[i][last] == vec[i+1][0]){
                vec[i].pop_back();
            }
        }
    }

template <digest::BadCharPolicy P, class T>
void thread_wind(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_wind<P, T>(thread_count, vec, seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h);
    }

template <digest::BadCharPolicy P, class T>
void thread_wind(unsigned thread_count, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& vec, 
    const std::string& seq, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_wind<P, T>(thread_count, vec, seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h);
    }

template <digest::BadCharPolicy P, class T>
std::vector<uint32_t> thread_wind_roll1(const char* seq, 
    size_t ind, unsigned k, uint32_t large_wind_kmer_am,
    digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am){
		std::vector<uint32_t> out;
        digest::WindowMin<P, T> dig(seq, ind + assigned_lwind_am + k + large_wind_kmer_am -1 -1, k, large_wind_kmer_am, ind, minimized_h);
        dig.roll_minimizer(assigned_lwind_am, out);
		return out;
    }

template <digest::BadCharPolicy P, class T>
std::vector<std::pair<uint32_t, uint32_t>> thread_wind_roll2(const char* seq, 
    size_t ind, unsigned k, uint32_t large_wind_kmer_am,
    digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am){
		std::vector<std::pair<uint32_t, uint32_t>> out;
        digest::WindowMin<P, T> dig(seq, ind + assigned_lwind_am + k + large_wind_kmer_am -1 -1, k, large_wind_kmer_am, ind, minimized_h);
        dig.roll_minimizer(assigned_lwind_am, out);
		return out;
    }

template <digest::BadCharPolicy P, class T>
void thread_sync(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_lwinds = (int)len - (int)start - (int)(k+large_wind_kmer_am)+2;
        if(large_wind_kmer_am == 0 || k < 4 ||start >= len || num_lwinds < 0 || (unsigned)num_lwinds < thread_count){
            throw BadThreadOutParams();
        }
        unsigned lwinds_per_thread = num_lwinds/thread_count;
        unsigned extras = num_lwinds % thread_count;
        vec.reserve(thread_count);
        std::vector<std::future<std::vector<uint32_t>>> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_lwind_am = lwinds_per_thread;
            if(extras > 0){
                ++(assigned_lwind_am);
                extras--;
            }

            thread_vector.emplace_back(std::async(thread_sync_roll1<P, T>,
                seq, ind, k, large_wind_kmer_am, minimized_h, assigned_lwind_am));

            ind += assigned_lwind_am;
        }
        for(auto& t: thread_vector)
        {

						vec.emplace_back(t.get());
        }
        
    }
template <digest::BadCharPolicy P, class T>
void thread_sync(unsigned thread_count, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& vec, 
    const char* seq, size_t len, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        int num_lwinds = (int)len - (int)start - (int)(k+large_wind_kmer_am)+2;
        if(large_wind_kmer_am == 0 || k < 4 ||start >= len || num_lwinds < 0 || (unsigned)num_lwinds < thread_count){
            throw BadThreadOutParams();
        }
        unsigned lwinds_per_thread = num_lwinds/thread_count;
        unsigned extras = num_lwinds % thread_count;
        vec.reserve(thread_count);
        std::vector<std::future<std::vector<std::pair<uint32_t, uint32_t>>>> thread_vector;

        size_t ind = start;
        for(unsigned i = 0; i < thread_count; i++){
            // issue is here
            // this will lead to a leak
            unsigned assigned_lwind_am = lwinds_per_thread;
            if(extras > 0){
                ++(assigned_lwind_am);
                extras--;
            }

            thread_vector.emplace_back(std::async(thread_sync_roll2<P, T>,
                seq, ind, k, large_wind_kmer_am, minimized_h, assigned_lwind_am));

            ind += assigned_lwind_am;
        }
        for(auto& t: thread_vector)
        {
			vec.emplace_back(t.get());
        }
        
    }

template <digest::BadCharPolicy P, class T>
void thread_sync(unsigned thread_count, std::vector<std::vector<uint32_t>>& vec, 
    const std::string& seq, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_sync<P, T>(thread_count, vec, seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h);
    }

template <digest::BadCharPolicy P, class T>
void thread_sync(unsigned thread_count, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& vec, 
    const std::string& seq, unsigned k, uint32_t large_wind_kmer_am, size_t start, 
    digest::MinimizedHashType minimized_h){
        thread_sync<P, T>(thread_count, vec, seq.c_str(), seq.size(), k, large_wind_kmer_am, start, minimized_h);
    }
    
template <digest::BadCharPolicy P, class T>
std::vector<uint32_t> thread_sync_roll1(const char* seq, 
    size_t ind, unsigned k, uint32_t large_wind_kmer_am,
    digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am){
		std::vector<uint32_t> out;
        digest::Syncmer<P, T> dig(seq, ind + assigned_lwind_am + k + large_wind_kmer_am -1 -1, k, large_wind_kmer_am, ind, minimized_h);
        dig.roll_minimizer(assigned_lwind_am, out);
		return out;
    }

template <digest::BadCharPolicy P, class T>
std::vector<std::pair<uint32_t, uint32_t>> thread_sync_roll2(const char* seq, 
    size_t ind, unsigned k, uint32_t large_wind_kmer_am,
    digest::MinimizedHashType minimized_h, unsigned assigned_lwind_am){
		std::vector<std::pair<uint32_t, uint32_t>> out;
        digest::Syncmer<P, T> dig(seq, ind + assigned_lwind_am + k + large_wind_kmer_am -1 -1, k, large_wind_kmer_am, ind, minimized_h);
        dig.roll_minimizer(assigned_lwind_am, out);
		return out;
    }
}
