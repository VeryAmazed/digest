#include <digest/window_minimizer.hpp>
#include <digest/syncmer.hpp>
#include <digest/mod_minimizer.hpp>
#include <variant>

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>> window_minimizer (
        const std::string &seq, unsigned k, unsigned large_window, bool include_hash=false) {
    digest::WindowMin<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive> digester (seq, k, large_window);
    if (include_hash) {
        std::vector<std::pair<uint32_t, uint32_t>> output;
        digester.roll_minimizer(seq.length(), output);
        return output;
    }
    else {
        std::vector<uint32_t> output;
        digester.roll_minimizer(seq.length(), output);
        return output;
    }
}
//std::vector<std::pair<size_t, size_t>> output;

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>> modimizer (
        const std::string &seq, unsigned k, uint32_t mod, bool include_hash=false) {
    digest::ModMin<digest::BadCharPolicy::SKIPOVER> digester (seq, k, mod);
    if (include_hash) {
        std::vector<std::pair<uint32_t, uint32_t>> output;
        digester.roll_minimizer(seq.length(), output);
        return output;
    }
    else {
        std::vector<uint32_t> output;
        digester.roll_minimizer(seq.length(), output);
        return output;
    }
}

std::variant<std::vector<uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>> syncmer (
        const std::string &seq, unsigned k, unsigned large_window, bool include_hash=false) {
    digest::Syncmer<digest::BadCharPolicy::SKIPOVER, digest::ds::Adaptive> digester (seq, k, large_window);
    if (include_hash) {
        std::vector<std::pair<uint32_t, uint32_t>> output;
        digester.roll_minimizer(seq.length(), output);
        return output;
    }
    else {
        std::vector<uint32_t> output;
        digester.roll_minimizer(seq.length(), output);
        return output;
    }
}