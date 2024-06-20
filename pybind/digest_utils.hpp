#include <digest/window_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <digest/mod_minimizer.hpp>

std::vector<uint32_t> window_minimizer(const std::string &seq, unsigned k, unsigned large_window) {
    digest::WindowMin<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive> digester (seq, k, large_window);
    std::vector<uint32_t> output;
    digester.roll_minimizer(seq.size(), output);
    return output;
}

std::vector<uint32_t> modimizer(const std::string &seq, unsigned k, uint32_t mod) {
    digest::ModMin<digest::BadCharPolicy::SKIPOVER> digester (seq, k, mod);
    std::vector<uint32_t> output;
    digester.roll_minimizer(seq.size(), output);
    return output;
}

std::vector<uint32_t> syncmer(const std::string &seq, unsigned k, unsigned large_window) {
    digest::Syncmer<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive> digester (seq, k, large_window);
    std::vector<uint32_t> output;
    digester.roll_minimizer(seq.size(), output);
    return output;
}