# Digest
C++ library which supports various minimizer schemes for digestion of DNA sequences  

# Implementation
Supports Mod Minimizers, Window Minimizers, and Syncmers  

Uses the cyclic or hash provided by [ntHash](https://github.com/bcgsc/ntHash). 
Tests were written using the [Catch2](https://github.com/catchorg/Catch2) unit testing framework.    
Benchmarking is done using Google's [benchmark](https://github.com/google/benchmark) library.

Non-ACTG character's cannot be hashed and must be handled using one of the two bad character handling policies. More details in the documentation.

Mod Minimzer classifies a kmer as a minimizer if the hash of the kmer is congruent to the user specified value in the user specified mod-space.  

Window Minimizer classifies a kmer as a minimizer if it is the smallest in the user specifed large window, using rightmost kmer to break ties.  

Syncmer classifies a large window as a minimizer if its smallest value is equal to the value of the hashes of the leftmost or rightmost kmer in the window (doesn't care if the smallest hash value is not unique). Note that because of how the large window is defined if you are using the SKIPOVER policy and your sequence has non-ACTG characters, it is possible for this large window to have varying lengths in terms of number of characters.  

# Install
We use [Meson](https://mesonbuild.com). (Very) old version will not work.

PREFIX is an absolute path to install location. If excluded, will install to system libraries.
```bash
meson setup --prefix=PREFIX --buildtype=release build
meson install -C build
```
This will generate `include` and `lib` folders.

# Usage
[Documentation](https://veryamazed.github.io/digest/)
* Digest objects require that the input string is kept in memory, unmodified.
* requires `c++17`

**Note: Do not use the new_seq functions, just create new objects. They currently don't work correctly.**

# Example
```cpp
digest::WindowMin<digest::ds::Naive<8>> wm(str, 16, 8);
std::vector<size_t> temp;
wm.roll_minimizer(100000, temp);
```
Example snippet to collect up to 100000 indices of minimizers.
A vector must be passed in, which will be appended to.
Each WindowMin / Syncmer object is templated by the algorithm / data structure to find minimizers.

A complete example and cli can be found [here](https://github.com/BenLangmead/gester/tree/main)

# Contributing
Use clang format version 17.  
run `ninja clang-format` before submitting a PR.

# Benchmark / Tests
```bash
meson setup build
cd build && meson compile
```
this will generate proper executables for benchmark/testing
