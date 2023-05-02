cmake_minimum_required(VERSION 3.5)

project(digester_test)
set (CMAKE_CXX_STANDARD 11)
############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(ntHash STATIC 
    src/ntHash/nthash.cpp src/ntHash/nthash_consts.cpp src/ntHash/nthash_lowlevel.cpp src/ntHash/status.cpp
)

target_include_directories(ntHash
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include/ntHash
)

add_library(Digester STATIC 
    src/Digester/mod_minimizer.cpp src/Digester/digester.cpp src/Digester/window_minimizer.cpp
)

target_include_directories(Digester
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include/Digester ${PROJECT_SOURCE_DIR}/include/ntHash
)

############################################################
# Catch2 Stuff
############################################################

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

############################################################
# Create an executable
############################################################

# Add an executable with the above sources

add_executable(test 
    test.cpp
)

# link the new hello_library target with the hello_binary target
target_link_libraries(Digester 
    PRIVATE
        ntHash
)

target_link_libraries( test
    PRIVATE 
        ntHash
    PRIVATE
        Digester
    PRIVATE
        Catch2::Catch2WithMain
)