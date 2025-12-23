#!/usr/bin/env bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export CMAKE_PREFIX_PATH="/opt/homebrew/opt/openblas"
cmake -S . -B build && cmake --build build -j
