#!/usr/bin/env bash
if [[ "$1" == [Rr] ]]; then
    mode=R
elif [[ "$1" == [Dd] ]]; then
    mode=D
elif [[ "$1" == [Tt] ]]; then
    mode=T
else
  read -p "mode ([R]elease/[D]ebug/[T]est): " mode
fi

if [[ $mode == [Rr] ]]; then
    mode=Release
    echo Release mode
elif [[ $mode == [Dd] ]]; then
    mode=Debug
    echo Debug mode
elif [[ $mode == [Tt] ]]; then
    mode=Test
    echo Test mode
else
    echo Try agin...
    exit
fi

rm -rf build
mkdir build
cd build

if [[ $mode == "Test" ]]; then
    cmake -DCMAKE_BUILD_TYPE=Debug -DTEST=on -S ..
else
  cmake -DCMAKE_BUILD_TYPE="$mode" -S ..
fi
make -j 8
cd ..
