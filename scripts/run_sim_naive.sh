#!/usr/bin/env bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"

cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DDTT_BUILD_EXAMPLES=ON >/dev/null
cmake --build "${BUILD_DIR}" --target dtt_sim_naive -j >/dev/null

N=${1:-512}
STEPS=${2:-20}
DT=${3:-0.01}
CUTOFF=${4:-"nan"}
GRAVITY=${5:-9.0}
DIST=${6:-"uniform"}
OUT_DIR="${ROOT_DIR}/output_naive"

# Clean existing output to avoid mixed legacy (.vtk) and XML (.vtp) frames.
if [[ -d "${OUT_DIR}" ]]; then
  rm -rf "${OUT_DIR}"
fi

"${BUILD_DIR}/dtt_sim_naive" "${N}" "${STEPS}" "${DT}" "${DIST}" "${CUTOFF}" "${GRAVITY}" "${OUT_DIR}"

echo "Done. Load ${OUT_DIR}/frames.pvd in ParaView."
