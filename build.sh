#!/usr/bin/env bash
####################################################################################################
# Copyright (c) 2021-2023. All rights reserved.                                                    #
# @Author: Zhucong Xi                                                                              #
# @Date: 12/4/21 6:42 AM                                                                           #
# @Last Modified by: zhucongx                                                                      #
# @Last Modified time: 8/22/23 11:40 AM                                                            #
####################################################################################################

# Ask for the mode and compiler if they are not passed as arguments
mode=$1
[[ -z $mode ]] && {
    read -p "Enter the mode - [R]elease, [D]ebug: " mode
}
compiler=$2
[[ -z $compiler ]] && {
    read -p "Enter the compiler - [d]efault, [g]cc, [i]ntel, [o]neapi, [c]lang, [m]pi: " compiler
}
flag=$3

mode=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
compiler=$(echo "$compiler" | tr '[:upper:]' '[:lower:]')

# Decide build mode based on user input
case $mode in
  r|release) mode=Release ;;
  d|debug) mode=Debug ;;
  *) echo "Invalid mode. Try again..."; exit 1 ;;
esac
# Decide build compiler based on user input
case $compiler in
  d|default) compiler_C=cc compiler_CXX=c++ ;;
  g|gcc) compiler_C=gcc compiler_CXX=g++ ;;
  i|intel) compiler_C=icc compiler_CXX=icpc ;;
  o|intel-oneapi) compiler_C=icx compiler_CXX=icpx ;;
  c|clang) compiler_C=clang compiler_CXX=clang++ ;;
  m|mpi) compiler_C=mpicc compiler_CXX=mpicxx ;;
  *) echo "Invalid compiler. Try again..." ; exit 1 ;;
esac

echo "Building in $mode mode using $compiler_C, $compiler_CXX compiler and additional flag ($flag)."
if [ -d "cmake-build" ]; then
  rm -rf cmake-build
fi
mkdir cmake-build; cd cmake-build || {
  echo "Could not create/access build directory. Check if you have the required permissions and try again."
  exit 1
}

# Pass the selected parameters to cmake
cmake -D CMAKE_C_COMPILER="$compiler_C" -D CMAKE_CXX_COMPILER="$compiler_CXX" \
-D CMAKE_BUILD_TYPE="$mode" -D ADDITIONAL_FLAGS="$flag"  -S ..
cmake --build . -j 8

# Moving back to the original directory
cd ..