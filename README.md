# ✂️ Digest: fast, multi-use $k$-mer sub-sampling library

<p align="left">
  <img width="900" alt="image1" src="https://github.com/oma219/digest/assets/32006908/09523db6-fd0b-49de-8e57-0d2cedef2a26">
  <br>
  <em>Visualization of different minimizer schemes supported in Digest and code example using library </em>
</p>


## What is the Digest library?
- a `C++` library that supports various sub-sampling schemes for $k$-mers in DNA sequences.
    - `Digest` library utilizes the rolling hash-function from [ntHash](https://github.com/bcgsc/ntHash) to order the $k$-mers in a window.
- a set of Python bindings that allow the user to run functions from the C++ library in Python.
  
## How to install and build into your project?
<img width="600" alt="image2" src="https://github.com/oma219/digest/assets/32006908/7cea427e-c22a-4271-a234-a2aafeb45413">

### Option 1: conda installation
Digest is available on bioconda. This installs both the C++ library and python library. The `include` and `lib` directories are in the conda environment dir (you can find it using `conda env list`).
```bash
conda install -c bioconda digest
```

### Option 2: Install from source
After cloning from GitHub, we use the [Meson](https://mesonbuild.com) build-system to install the library. 
- `PREFIX` is an absolute path to library files will be install (`*.h` and `*.a` files)
    - **IMPORTANT**: `PREFIX` should not be the root directory of the `digest/` repo to avoid any issues with installation.
    - We suggest using `--prefix=$(pwd)/build` from within the root directory of the `digest/` repo.
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

To use Digest in your C++ project, you just need to include the header files (`*.h`) and library file (`*.a`) that were installed in the first step. Assuming that `build/` is the directory you installed them in, here is how you can compile.

```bash
g++ -std=c++17  -o main main.cpp -I build/include/ -L build/lib -lnthash
```

## Detailed Look at Example Usage (2 ways):

There are three types of minimizer schemes that can be used:

1. **Windowed Minimizer**: classifies a kmer as a minimizer if it is the smallest in the user specifed large window, using rightmost kmer to break ties.
2. **Modimizer**: classifies a kmer as a minimizer if the hash of the kmer is congruent to the user specified value in the user specified mod-space.
3. **Syncmer**: classifies a large window as a minimizer if its smallest value is equal to the value of the hashes of the leftmost or rightmost kmer in the window (doesn't care if the smallest hash value is not unique). Note that because of how the large window is defined if you are using the SKIPOVER policy and your sequence has non-ACTG characters, it is possible for this large window to have varying lengths in terms of number of characters.

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

### Documentation:
Documentation generated with Doxygen can be found [here](https://veryamazed.github.io/digest/)

## Python binding support

Included in the library are function bindings for each sub-sampling scheme for use in Python. The simplest way to install the python module is through conda (`conda install bioconda::digest`). To install the Python module from source, first install the library with `meson` (see above for detailed instructions), and install with `pip`. For this setup, the `meson` prefix must be set to `--prefix=/$DIGEST_REPO/build`:
```
meson setup --prefix=$(pwd)/build --buildtype=release build
meson install -C build
pip install .
```
Alternatively, copy the `lib` and `include` directories from the earlier meson installation to a directory in the repo called `build`, and run `pip install .`

We recommend using a conda or python virtual environment.
Once installed, you can import and use the Digest library in Python:
```
>>> from digest import window_minimizer, syncmer, modimizer
>>> window_minimizer('ACGTACGTAGCTAGCTAGCTAGCTGATTACATACTGTATGCAAGCTAGCTGATCGATCGTAGCTAGTGATGCTAGCTAC', k=5, w=11)
[4, 5, 16, 19, 21, 26, 27, 35, 39, 49, 57, 63, 68]
>>> modimizer('ACGTACGTAGCTAGCTAGCTAGCTGATTACATACTGTATGCAAGCTAGCTGATCGATCGTAGCTAGTGATGCTAGCTAC', k=5, mod=5)
[23, 34, 38, 40, 62, 67]
>>> syncmer('ACGTACGTAGCTAGCTAGCTAGCTGATTACATACTGTATGCAAGCTAGCTGATCGATCGTAGCTAGTGATGCTAGCTAC', k=5, w=15)
[0, 3, 4, 5, 7, 12, 13, 27, 35, 49]
>>> modimizer('ATCGTGCATCA', k=4, mod=2, include_hash=True)
[(0, 1122099596), (2, 249346952), (4, 227670418), (7, 123749036)]
>>> seq = 'ACGTACGTAGCTAGCTAGCTAGCTGATTACATACTGTATGCAAGCTAGCTGATCGATCGTAGCTAGTGATGCTAGCTAC'
>>> [seq[p:p+5] for p in window_minimizer(seq, k=5, w=11)]
['ACGTA', 'CGTAG', 'AGCTA', 'TAGCT', 'GCTGA', 'TTACA', 'TACAT', 'GTATG', 'GCAAG', 'TGATC', 'CGTAG', 'TAGTG', 'ATGCT']
```
## CLI
A basic cli can be found [here](https://github.com/BenLangmead/gester/tree/main)

## Benchmark / Tests
```bash
cd build
meson configure -Dbuildtype=debug
meson compile
```
This will set the build flage from release to debug allowing you to generate proper executables for benchmark/testing. The executables will be located in the build folder and can be run directly from there. You can look at the meson.build file for more details.

## Contributing
Use clang format version 17.  
run `ninja clang-format` before submitting a PR.

We have also implemented parallel execution in the python library:
```
>>> import timeit
>>> timeit.timeit('window_minimizer(seq, k=5, w=11)', 
...     setup='from digest import window_minimizer; seq = open("tests/bench/chrY.fa", "rb").read()',
...     number=10)
19.85955114942044
>>> timeit.timeit('window_minimizer(seq, k=5, w=11, num_threads=4)',
...     setup='from digest import window_minimizer; seq = open("tests/bench/chrY.fa", "rb").read()',
...     number=10)
10.327348310500383
```

Functions can also return numpy arrays with the `*_np` function variant:
```
>>> from digest import window_minimizer_np
>>> window_minimizer_np(b'ACGTACGTAGCTAGCTAGCTAGCTGATTACATACTGTATGCAAGCTAGCTGATCGATCGTAGCTAGTGATGCTAGCTAC', k=5, w=11)
array([ 4,  5, 16, 19, 21, 26, 27, 35, 39, 49, 57, 63, 68], dtype=uint32)
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
