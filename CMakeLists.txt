####################################################################################################
# Copyright (c) 2020-2023. All rights reserved.                                                    #
# @Author: Zhucong Xi                                                                              #
# @Date: 1/16/20 1:55 AM                                                                           #
# @Last Modified by: zhucongx                                                                      #
# @Last Modified time: 8/21/23 3:33 PM                                                             #
####################################################################################################

cmake_minimum_required(VERSION 3.14)
#################################################### SET  COMPILERS ####################################################
#set(CMAKE_C_COMPILER "icc") # {gcc, gcc-13, mpicc, icc, mpiicc, clang-16}
#set(CMAKE_CXX_COMPILER "icpc") # {g++, g++-13, mpicxx, icpc, mpiicpc, clang-16}
#set(CMAKE_Fortran_COMPILER "ifort") # {gfortran, gfortran-13, mpifort, ifort, mpiifort, flang}
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)

project(lmc VERSION 0.1.0)
set(CMAKE_VERBOSE_MAKEFILE ON)
###################################################### SET  PATHS ######################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
###################################################### SET  FLAGS ######################################################
include(CheckCXXCompilerFlag)
macro(add_cxx_compiler_flags FLAGS)
    separate_arguments(FLAGS_LIST UNIX_COMMAND "${FLAGS}")
    foreach (FLAG ${FLAGS_LIST})
        string(REGEX REPLACE "-" "" SFLAG1 ${FLAG})
        string(REGEX REPLACE "\\+" "p" SFLAG ${SFLAG1})
        check_cxx_compiler_flag(${FLAG} COMPILER_SUPPORT_${SFLAG})
        if (COMPILER_SUPPORT_${SFLAG})
            if (CMAKE_CXX_FLAGS STREQUAL "")
                set(CMAKE_CXX_FLAGS "${FLAG}")
            else ()
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
            endif ()
        endif ()
    endforeach ()
endmacro()
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_cxx_compiler_flags("-g -O0")
    # add_cxx_compiler_flags("-Werror")
    add_cxx_compiler_flags("-pedantic -Wall -Wextra -Wundef -Wcast-align -Wcast-qual -Wctor-dtor-privacy")
    add_cxx_compiler_flags("-Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wextra-semi")
    add_cxx_compiler_flags("-Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast")
    add_cxx_compiler_flags("-Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo")
    add_cxx_compiler_flags("-Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wno-unused -Weffc++")
    add_cxx_compiler_flags("-Wchar-subscripts -Wnon-virtual-dtor -Wunused-local-typedefs -Wpointer-arith")
    add_cxx_compiler_flags("-Wwrite-strings -Wformat-security -Wshorten-64-to-32 -Wenum-conversion")
    add_cxx_compiler_flags("-Wc++11-extensions -Wdouble-promotion -Wconversion -Wunused -Winline -Wunreachable-code")
    add_cxx_compiler_flags("-Wmisleading-indentation -pedantic-errors -Wnull-dereference -Wfatal-errors")
    add_cxx_compiler_flags("-wd593 -wd981 -wd1599 -wd2012 -wd2014 -wd2015 -wd2017 -wd2021 -wd2022 -wd2027")
    add_cxx_compiler_flags("-wd2047 -wd2048 -wd2304 -wd2338 -wd10441 -wd11071")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
    add_cxx_compiler_flags("-O3")
    add_cxx_compiler_flags("-DNDEBUG -DEIGEN_NO_DEBUG")
    add_cxx_compiler_flags("-fno-math-errno -fno-trapping-math -ffp-contract=fast")
#    add_cxx_compiler_flags("-ffast-math") # NOTE: this flag may lead to precision loss
#    add_cxx_compiler_flags("-march=native") # NOTE: this flag better to be added using ADDITIONAL_FLAGS
endif ()
if (DEFINED ADDITIONAL_FLAGS AND NOT "${ADDITIONAL_FLAGS}" STREQUAL "")
    add_cxx_compiler_flags("${ADDITIONAL_FLAGS}")
endif ()
message(STATUS "Detecting ${CMAKE_CXX_COMPILER_ID} compiler")
message(STATUS "Current build type is ${CMAKE_BUILD_TYPE}")
###################################################### BUILD SRCS ######################################################
file(GLOB PROJECT_SRC ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*/src/*.cpp)
#file(GLOB PROJECT_SRC
#     ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/cfg/src/*.cpp
#     ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/pred/src/*.cpp
#     # ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/ce/src/*.cpp
#     ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/mc/src/*.cpp
#     ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/ansys/src/*.cpp
#     ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/api/src/*.cpp
#)

file(GLOB PROJECT_SOURCE ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*/src)
foreach (dir ${PROJECT_SOURCE})
    message(STATUS "Found source directory: ${dir}")
endforeach ()
file(GLOB PROJECT_INCLUDE ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/*/include)
foreach (dir ${PROJECT_INCLUDE})
    message(STATUS "Found include directory: ${dir}")
endforeach ()

add_executable(${PROJECT_NAME} ${PROJECT_NAME}/main.cpp ${PROJECT_SRC})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_INCLUDE})
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}.exe")
###################################################### BUILD LIBS ######################################################
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW) # Use modern timestamp behavior for FetchContent
    set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)
endif()
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW) # Use modern FindBoost module
    set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
endif ()
if (POLICY CMP0146)
    cmake_policy(SET CMP0146 NEW) # Use modern FindCUDA module
    set(CMAKE_POLICY_DEFAULT_CMP0146 NEW)
endif()
include(FetchContent)
function(target_link_libraries_system target scope)
    set(libs ${ARGN})
    foreach (lib ${libs})
        get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
        if (lib_include_dirs AND NOT lib_include_dirs STREQUAL "NOTFOUND")
            target_include_directories(${target} SYSTEM ${scope} ${lib_include_dirs})
        endif ()
        target_link_libraries(${target} ${scope} ${lib})
    endforeach ()
endfunction()
#function(declare_system_library target)
#    message(STATUS "Declaring system library ${target}")
#    get_target_property(target_aliased_name ${target} ALIASED_TARGET)
#    if (target_aliased_name)
#        set(target ${target_aliased_name})
#    endif ()
#    set_target_properties(${target} PROPERTIES
#                          INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
#                          $<TARGET_PROPERTY:${target},INTERFACE_INCLUDE_DIRECTORIES>)
#endfunction()
####################################################### PKGCONFIG ######################################################
#find_package(PkgConfig)
######################################################## BOOST #########################################################
#set(Boost_USE_DEBUG_LIBS OFF) # ignore debug libs and
#set(Boost_USE_RELEASE_LIBS ON) # only find release libs
find_package(Boost REQUIRED COMPONENTS filesystem iostreams)
target_link_libraries_system(${PROJECT_NAME} PRIVATE Boost::boost Boost::filesystem Boost::iostreams)
######################################################### MPI ##########################################################
find_package(MPI REQUIRED)
target_link_libraries_system(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
######################################################## OPENMP ########################################################
find_package(OpenMP REQUIRED)
target_link_libraries_system(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
######################################################## EIGEN #########################################################
find_package(MKL)
if (MKL_FOUND AND NOT (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    message(STATUS "MKL found. Using MKL BLAS/LAPACK.")
    set(BLA_VENDOR Intel10_64ilp)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    target_link_libraries_system(${PROJECT_NAME} PRIVATE MKL::MKL)
else ()
    message(STATUS "MKL not found.")
    if (DEFINED ENV{MKLROOT})
        message(STATUS "But MKLROOT found, try MKL again: $ENV{MKLROOT}")
        set(BLA_VENDOR Intel10_64ilp)
        add_definitions(-DEIGEN_USE_MKL_ALL)
        if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
            add_definitions(-mkl)
        endif ()
    else ()
        message(STATUS "MKLROOT not found. Using generic BLAS/LAPACK.")
        set(BLA_VENDOR Generic)
        add_definitions(-DEIGEN_USE_LAPACK -DEIGEN_USE_BLAS)
    endif ()
    set(BLA_VENDOR Generic)
    find_package(BLAS)
    find_package(LAPACK)
    if (BLAS_FOUND AND LAPACK_FOUND)
        target_link_libraries_system(${PROJECT_NAME} PRIVATE BLAS::BLAS)
        target_link_libraries_system(${PROJECT_NAME} PRIVATE LAPACK::LAPACK)
    else ()
        message(STATUS "BLAS/LAPACK not found. Using Eigen3 only.")
    endif ()
endif ()

find_package(Eigen3)
# If a suitable system version wasn't found, download the exact version we need.
if (NOT Eigen3_FOUND AND NOT EIGEN3_FOUND)
    message(STATUS "Eigen3 not found.")
    message(STATUS "Downloading Eigen3 from Githab")
    FetchContent_Declare(
            Eigen3
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 5.0.0
            # GIT_PROGRESS TRUE
            # GIT_SHALLOW TRUE
    )
#    set(FETCHCONTENT_QUIET ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(Eigen3)
elseif (Eigen3_FOUND)
    message(STATUS "Eigen3 found. Using version ${Eigen3_VERSION}")
elseif (EIGEN3_FOUND)
        message(STATUS "Eigen3 found. Using version ${EIGEN3_VERSION}")
endif ()
# Link to the Eigen3::Eigen target, which is now available either from the system or FetchContent.
target_link_libraries_system(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
###################################################### ARMADILLO #######################################################
#find_package(ARMADILLO REQUIRED)
#target_include_directories(${PROJECT_NAME} PRIVATE
#        ${ARMADILLO_INCLUDE_DIRS})
#target_link_libraries_system(${PROJECT_NAME} PRIVATE
#        ${ARMADILLO_LIBRARIES})
###################################################### TENSORFLOW ######################################################
#find_package(TensorFlow REQUIRED)
#target_include_directories(${PROJECT_NAME} PRIVATE
#        ${TensorFlow_INCLUDE_DIRS})
#target_link_libraries_system(${PROJECT_NAME} PRIVATE
#        ${TensorFlow_LIBRARIES})
######################################################### JSON #########################################################
# Typically you don't care so much for a third party library's tests to be run from your own project's code.
# set(JSON_BuildTests OFF CACHE INTERNAL "")
# If you only include this third party in PRIVATE source files, you do not
# need to install it when your main project gets installed.
# set(JSON_Install OFF CACHE INTERNAL "")
# add_subdirectory(external/json)
message(STATUS "Downloading nlohmann_json from Github")
FetchContent_Declare(
        json
        URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
)
FetchContent_MakeAvailable(json)
target_link_libraries_system(${PROJECT_NAME} PRIVATE nlohmann_json::nlohmann_json)
######################################################## SPGLIB ########################################################
# message(STATUS "Downloading spglib from Github")
# FetchContent_Declare(
#         Spglib
#         GIT_REPOSITORY https://github.com/spglib/spglib.git
#         GIT_TAG v2.1.0-rc2
#         GIT_PROGRESS TRUE
#         GIT_SHALLOW TRUE
#         SYSTEM TRUE
# )

# set(CMAKE_BUILD_TYPE_BACKUP ${MAKE_BUILD_TYPE})
# set(CMAKE_BUILD_TYPE Release)
# set(SPGLIB_USE_OMP ON)
# set(SPGLIB_WITH_TESTS OFF)
# FetchContent_MakeAvailable(Spglib)
# set(CMAKE_BUILD_TYPE ${MAKE_BUILD_TYPE_BACKUP})

# #FetchContent_GetProperties(Spglib)
# #if (NOT Spglib_POPULATED)
# #    FetchContent_Populate(Spglib)
# #    add_subdirectory(${spglib_SOURCE_DIR} ${spglib_BINARY_DIR} EXCLUDE_FROM_ALL)
# #endif ()
# target_link_libraries_system(${PROJECT_NAME} PRIVATE Spglib::symspg)
###################################################### GOOGLETEST ######################################################
#set(TEST "on")
#IF(TEST)
#    message(STATUS "TEST mode on")
#    include(FetchContent)
#    FetchContent_Declare(
#            googletest
#            GIT_REPOSITORY https://github.com/google/googletest/
#            GIT_TAG release-1.11.0)
#    # For Windows: Prevent overriding the parent project's compiler/linker settings
#    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#    FetchContent_MakeAvailable(googletest)
#
#    enable_testing()
#    include(GoogleTest)
##    macro(package_add_test TESTNAME)
##        add_executable(${TESTNAME} ${ARGN})
##        target_link_libraries(${TESTNAME} gtest gmock gtest_main)
##        gtest_discover_tests(${TESTNAME}
##                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
##                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test"
##                )
##        set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
##    endmacro()
#
##    macro(package_add_test_with_libraries TESTNAME FILES LIBRARIES)
##        add_executable(${TESTNAME} ${FILES})
##        target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})
##        gtest_discover_tests(${TESTNAME}
##                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
##                PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/test"
##                )
##        set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
##    endmacro()
#    set(TEST_NAME my_test)
#    add_executable(${TEST_NAME} gtest.cpp ${PROJECT_SRC})
#    set_target_properties(${TEST_NAME} PROPERTIES
#            OUTPUT_NAME "${TEST_NAME}.exe")
#    target_include_directories(${TEST_NAME} PRIVATE
#            ${PROJECT_INCLUDE}
#            ${Boost_INCLUDE_DIR})
#
#    target_link_libraries(${TEST_NAME} gtest gmock gtest_main )
#    gtest_discover_tests(${TEST_NAME})
#ELSE()
#    message(STATUS "TEST mode off")
#ENDIF()
###################################################### SANITIZER #######################################################
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
#add_cxx_compiler_flags("-fsanitize=address -fno-omit-frame-pointer") # address
#add_cxx_compiler_flags("-fsanitize=memory -fPIE -pie -g") # memory
#add_cxx_compiler_flags("-fsanitize=leak") # leak
#add_cxx_compiler_flags("-fsanitize=thread -fPIE -pie -g") # thread
#add_cxx_compiler_flags("-fsanitize=undefined") # undefined
##################################################### SHOW SUMMARY #####################################################
MESSAGE(STATUS "====================================================================================================")
MESSAGE(STATUS "Project Name:     ${PROJECT_NAME}")
MESSAGE(STATUS "C Compiler:       ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
MESSAGE(STATUS "C Standard:       C ${CMAKE_C_STANDARD}")
MESSAGE(STATUS "C++ Compiler:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
MESSAGE(STATUS "C++ Standard:     C++ ${CMAKE_CXX_STANDARD}")
MESSAGE(STATUS "Build type:       ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Build type flags: ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Output Path:      ${EXECUTABLE_OUTPUT_DIRECTORY}")
MESSAGE(STATUS "====================================================================================================")
