# Digest
C++ library which supports various minimizer schemes for digestion of DNA sequences  

# Implementation (Most of the documentation is in the code)
Supports Mod Minimizers, Window Minimizers, and Syncmers  

Uses the cyclic or hash provided by [ntHash](https://github.com/bcgsc/ntHash). For now I just downloaded the essential files off their github and compiled it myself but I may change how I link in ntHash in the future.  

Tests were written using the [Catch2](https://github.com/catchorg/Catch2) unit testing framework.  

As ntHash only hashes ACTG characters, this library can also only hash ACTG characters. Any kmer that contains a non-ACTG character will be skipped over, for example is the kmer size is 4, and the sequence is ACCGNATTGC, only ACCG, ATTG, and TTGC will be hashed and considered for being a minimizer.  

Mod Minimzer classifies a kmer as a minimizer if the hash of the kmer is congruent to the user specified value in the user specified mod-space.  

Window Minimizer classifies a kmer as a minimizer if it is the smallest in the user specifed large window, using rightmost kmer to break ties.  

Syncmer classifies a large window as a minimizer if its smallest value is equal to the value of the hashes of the leftmost or rightmost kmer in the window (doesn't care if the smallest hash value is not unique). Because of how the large window is defined and how this library handles skipping over non-ACTG characters, if your sequence has non-ACTG characters, it is possible for this large window to have varying lengths in terms of number of characters.  

# Install
We use [Meson](https://mesonbuild.com). Older version will not work.

PREFIX is an absolute path to install location.
```
meson setup --prefix=PREFIX build
meson install -C build
```
This will generate include/ and lib/ folders.

# Usage
* Headers at `#include <digest/___.hpp>`
* Classes are in `digest` namespace

# Todo
License  
expected density?  
combine headers  
