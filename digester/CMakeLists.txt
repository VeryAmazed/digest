cmake_minimum_required(VERSION 3.5)

project(Digester)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS_RELEASE "-O2")
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
    src/Digester/mod_minimizer.cpp src/Digester/digester.cpp src/Digester/window_minimizer.cpp src/Digester/syncmer.cpp
)

target_include_directories(Digester
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include/Digester ${PROJECT_SOURCE_DIR}/include/ntHash
)

target_link_libraries( Digester
    PRIVATE
        ntHash  
)

