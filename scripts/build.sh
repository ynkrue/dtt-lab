#!/usr/bin/env bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CC=/opt/homebrew/bin/gcc-13
export CXX=/opt/homebrew/bin/g++-13
cmake -S . -B build -DCBLAS_INCLUDE_DIR=/opt/homebrew/opt/openblas/include \
                    -DCMAKE_C_COMPILER=$CC \
                    -DCMAKE_CXX_COMPILER=$CXX \
                    -DCMAKE_BUILD_TYPE=Release

cmake --build build -- -j 8
