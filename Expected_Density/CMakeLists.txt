cmake_minimum_required(VERSION 3.5)

project(Expected)

set (CMAKE_CXX_STANDARD 11)

Include(FetchContent)

FetchContent_Declare(
  Digest
  GIT_REPOSITORY https://github.com/VeryAmazed/Digest
  GIT_TAG        0ebc52f89c5126a6d85c71224cb06251d9d43418
)

FetchContent_MakeAvailable(Digest)

add_executable(expected
    ACTG.cpp
)

target_link_libraries( expected
    PRIVATE
        digester
)

