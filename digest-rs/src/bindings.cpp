#include <digest/window_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <digest/mod_minimizer.hpp>
#include <digest/data_structure.hpp>
#include <cstring>
#include <memory>
#include <vector>
#include <string>

extern "C" {
    using namespace digest;

    // Window minimizer wrapper functions
    size_t window_minimizer(const char* seq, size_t len, unsigned k, unsigned large_window, uint32_t* out) {
        try {
            std::string sequence(seq, len);
            WindowMin<BadCharPolicy::SKIPOVER, ds::Adaptive> digester(sequence, k, large_window);
            std::vector<uint32_t> output;
            digester.roll_minimizer(sequence.length(), output);
            if (!output.empty()) {
                std::memcpy(out, output.data(), output.size() * sizeof(uint32_t));
            }
            return output.size();
        } catch (...) {
            return 0;
        }
    }
    
    // Modimizer wrapper function
    size_t modimizer(const char* seq, size_t len, unsigned k, uint32_t mod_val, uint32_t* out) {
        try {
            std::string sequence(seq, len);
            ModMin<BadCharPolicy::SKIPOVER> digester(sequence, k, mod_val);
            std::vector<uint32_t> output;
            digester.roll_minimizer(sequence.length(), output);
            if (!output.empty()) {
                std::memcpy(out, output.data(), output.size() * sizeof(uint32_t));
            }
            return output.size();
        } catch (...) {
            return 0;
        }
    }

    // Syncmer wrapper function
    size_t syncmer(const char* seq, size_t len, unsigned k, unsigned large_window, uint32_t* out) {
        try {
            std::string sequence(seq, len);
            Syncmer<BadCharPolicy::SKIPOVER, ds::Adaptive> digester(sequence, k, large_window);
            std::vector<uint32_t> output;
            digester.roll_minimizer(sequence.length(), output);
            if (!output.empty()) {
                std::memcpy(out, output.data(), output.size() * sizeof(uint32_t));
            }
            return output.size();
        } catch (...) {
            return 0;
        }
    }
} 