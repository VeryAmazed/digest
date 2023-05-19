# Digest
C++ library which supports various minimizer schemes for digestion of DNA sequences <br>

The project is functional. I'm not sure how optimized it is.
I'm not done with this project but I am putting it on hold for the summer. I'll also make this README look better when I finish.
# Implementation (Most of the documentation is in the code)
Supports Mod Minimizers, Window Minimizers, and Syncmers <br>

Uses the cyclic or hash provided by [ntHash](https://github.com/bcgsc/ntHash). For now I just downloaded the essential files off their github and compiled it myself but I may change how I link in ntHash in the future. <br>

Tests were written usig the [Catch2](https://github.com/catchorg/Catch2) unit testing framework. <br>

As ntHash only hashes ACTG characters, this library can also only hash ACTG characters. Any kmer that contains a non-ACTG character will be skipped over, for example is the kmer size is 4, and the sequence is ACCGNATTGC, only ACCG, ATTG, and TTGC will be hashed and considered for being a minimizer. <br>

Mod Minimzer classifies a kmer as a minimizer if the hash of the kmer is congruent to the user specified value in the user specified mod-space. <br>

Window Minimizer classifies a kmer as a minimizer if it is the smallest in the user specifed large window, using rightmost kmer to break ties. <br>

Syncmer classifies a large window as a minimizer if its smallest value is equal to the value of the hashes of the leftmost or rightmost kmer in the window (doesn't care if the smallest hash value is not unique). Because of how the large window is defined and how this library handles skipping over non-ACTG characters, if your sequence has non-ACTG characters, it is possible for this large window to have varying lengths in terms of number of characters. <br>
# Usage
```
FetchContent_Declare(
  Digest
  GIT_REPOSITORY https://github.com/VeryAmazed/Digest
  GIT_TAG        v0.1.1
)
```
Supports C++ 11 or later. Just copy the above code into your cmake file and that will generate a subdirectory called digester which you can then link into your project. Below is an example.
```
target_link_libraries( exec
    PRIVATE
        digester
)
```
You can actually look in the Expected_Density folder to see how I linked in this library. 
