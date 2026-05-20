#!/bin/bash

cd ./Mutationpp

mkdir build && cd build

cmake \
    -DCMAKE_INSTALL_PREFIX:PATH=$(realpath ../install) \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_DOCUMENTATION=ON \
    ..

make -j4 install
