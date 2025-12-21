#!/usr/bin/env bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cmake -S . -B build > /dev/null
cmake --build build -j > /dev/null
ctest --test-dir build
