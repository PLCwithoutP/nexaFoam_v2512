#!/bin/bash

echo "Including from: $MPP_INCLUDE"
echo "Linking from:   $MPP_LIB"

g++ -std=c++11 \
    -I$MPP_INCLUDE \
    -I/usr/include/eigen3 \
    mutationppMappingTest.cpp \
    -L$MPP_LIB \
    -lmutation++ \
    -Wl,-rpath,$MPP_LIB \
    -o test_mpp

./test_mpp
