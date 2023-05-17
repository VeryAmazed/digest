cmake_minimum_required(VERSION 3.5)

project(digester_test)

add_subdirectory(digester)

#message(STATUS "Testing. Did we get here 1")

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    #message(STATUS "Testing. Did we get here 2")
    
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


    target_link_libraries( test
        PRIVATE
            digester
        PRIVATE
            Catch2::Catch2WithMain
    )
endif()

