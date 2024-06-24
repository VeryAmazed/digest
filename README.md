# ✂️ Digest: fast, multi-use $k$-mer sub-sampling library

<p align="left">
  <img width="900" alt="image1" src="https://github.com/oma219/digest/assets/32006908/09523db6-fd0b-49de-8e57-0d2cedef2a26">
  <br>
  <em>Visualization of different minimizer schemes supported in Digest and code example using library </em>
</p>


## What is the Digest library?
- a `C++` library that supports various sub-sampling schemes for $k$-mers in DNA sequences.
    - `Digest` library utilizes the rolling hash-function from [ntHash](https://github.com/bcgsc/ntHash) to order the $k$-mers in a window.

## How to install and build into your project?
<img width="600" alt="image2" src="https://github.com/oma219/digest/assets/32006908/7cea427e-c22a-4271-a234-a2aafeb45413">

### Step 1: Install library

After cloning from GitHub, we use the [Meson](https://mesonbuild.com) build-system to install the library. 
- `PREFIX` is an absolute path to library files will be install (`*.h` and `*.a` files)
    - **IMPORTANT**: `PREFIX` should not be the root directory of the `Digest/` repo to avoid any issues with installation.
- These commands generate an `include` and `lib` folders in `PREFIX` folder

```bash
git clone https://github.com/VeryAmazed/digest.git

meson setup --prefix=<PREFIX> --buildtype=release build
meson install -C build
```

### Step 2: Include Digest in your project

#### (a) Using `Meson`: 

If your coding project uses `Meson` to build the executable(s), you can include a file called `subprojects/digest.wrap` in your repository and let Meson install it for you.

#### (b) Using `g++`:

To use Digest in your C++ project, you just need to include the header files (`*.h`) and library file (`*.a`) that were installed in the first step. Assuming that `install/` is the directory you installed them in, here is how you can compile.

```bash
g++ -std=c++17  -o main main.cpp -I install/include/ -L install/lib -lnthash
```

## Detailed Look at Example Usage (2 ways):

There are three types of minimizer schemes that can be used:

1. Windowed Minimizer
2. Modimizer
3. Syncmer

The general steps to use Digest is as follows: (1) include the relevant header files, (2) declare the Digest object and (3) find the positions where the minimizers are present in the sequence.

### 1. Find positions of minimizers:
```cpp
#include "digest/digester.hpp"
#include "digest/window_minimizer.hpp"

digest::WindowMin<digest::BadCharPolicy::WRITEOVER, digest::ds::Adaptive> digester (dna, 15, 7);

std::vector<size_t> output;
digester.roll_minimizer(100, output);
```
- This code snippet will find up to 100 Windowed Minimizers and store their positions in the vector called `output`.
- `digest::BadCharPolicy::WRITEOVER` means that anytime the code encounters an non-`ACTG` character, it will replace it with an `A`.
    - `digest::BadCharPolicy::SKIPOVER` will skip any $k$-mers with non-`ACTG` characters
- `digest::ds::Adaptive` is our recommended data-structure for finding the minimum value in a window (see wiki for other options)

### 2. Find both positions and hash values of minimizers
If you would like to obtain both the positions and hash values for each minimizer, you can pass a vector of paired integers to do so.

```
std::vector<std::pair<size_t, size_t>> output;
digester.roll_minimizer(100, output);
```


<!---
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

-->
