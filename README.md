# dtt-lab

CPU dual-tree traversal (DTT) sandbox for N-body style problems, with GPU extension (cuda kernels wip).
The library implements a naive baselines, spatial trees, and layer prune/accept/refine rules.
Translating the control structure to GPUs is in the works. Output snapshots will be written as VTK files so you can
visualize and animate runs in ParaView.

## Prerequisites

- CMake ≥ 3.26
- C++20 compiler
- BLAS
- Optional: `clang-format`, `clang-tidy`

## Configure, build, test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

### Useful options

- `-DDTT_BUILD_TESTS=ON` (default): build GoogleTest suite.
- `-DDTT_BUILD_BENCHMARKS=ON` (default): build Google Benchmark harnesses.
- `-DDTT_BUILD_EXAMPLES=ON` (default): build example executables.
- `-DDTT_ENABLE_CUDA=ON`: enable CUDA language for GPU.
- `-DDTT_ENABLE_CLANG_TIDY=ON`: run clang-tidy if found.

### Formatting

If `clang-format` is available:

```bash
cmake --build build --target format
```

### Benchmarks

```bash
cmake -S . -B build -DDTT_BUILD_BENCHMARKS=ON
cmake --build build --target dtt-bench
./build/dtt-bench
```

### Simulation example + ParaView

```bash
cmake -S . -B build -DDTT_BUILD_EXAMPLES=ON
cmake --build build --target dtt-sim
./build/dtt-sim 512 20 0.01 output
```

Open `output/frames.pvd` in ParaView and press Play.

## Project layout

- `include/` public headers
- `src/` library sources
- `cuda/` public headers and source code for cude device code
- `tests/` GoogleTest suite
- `benchmarks/` performance benchmark harnesses
- `examples/` small runs that emit VTK output for ParaView
- `docs/CONVENTIONS.md` coding standards and design notes
