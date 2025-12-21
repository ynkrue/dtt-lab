#!/usr/bin/env bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"

if ! command -v clang-format >/dev/null 2>&1; then
  echo "clang-format not found in PATH. Please install it first." >&2
  exit 1
fi

cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" >/dev/null
cmake --build "${BUILD_DIR}" --target format
