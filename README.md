# Digest
C++ library which supports various minimizer schemes for digestion of DNA sequences

Not yet finished, but the code is documented and tested

Will be usable, but not fully optimized soon (hopefully) 

FetchContent_Declare(
  Digest
  GIT_REPOSITORY https://github.com/VeryAmazed/Digest
  GIT_TAG        0ebc52f89c5126a6d85c71224cb06251d9d43418
)

target_link_libraries( exec
    PRIVATE
        digester
)
