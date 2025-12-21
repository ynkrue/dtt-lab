# dtt-lab

CPU-first dual-tree traversal (DTT) sandbox for N-body style problems, with a clean path to CUDA.
We start with naive baselines, add spatial trees, and layer prune/accept/refine rules before
translating the control structure to GPUs. Output snapshots will be written as VTK files so you can
visualize and animate runs in ParaView.

## Prerequisites

- CMake ≥ 3.27
- C++20 compiler (Clang or GCC recommended)
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
- `-DDTT_ENABLE_CUDA=ON`: enable CUDA language for GPU (no kernels yet).
- `-DDTT_ENABLE_CLANG_TIDY=ON`: run clang-tidy if found.
- `-DDTT_BUILD_BENCHMARKS=ON` (default): build Google Benchmark harnesses.
- `-DDTT_BUILD_EXAMPLES=ON` (default): build example executables.

### Formatting

If `clang-format` is available:

```bash
cmake --build build --target format
```

### Benchmarks

```bash
cmake -S . -B build -DDTT_BUILD_BENCHMARKS=ON
cmake --build build --target dtt_bench
./build/dtt_bench
```

### Simulation example + ParaView

```bash
cmake -S . -B build -DDTT_BUILD_EXAMPLES=ON
cmake --build build --target dtt_sim
./build/dtt_sim 512 20 0.01 output
```
Open `output/frames.pvd` in ParaView and press Play.

## Project layout

- `include/` public headers (core data structures, traversal interfaces)
- `src/` library sources
- `tests/` GoogleTest suite
- `benchmarks/` performance harnesses
- `examples/` small runs that emit VTK output for ParaView
- `docs/CONVENTIONS.md` coding standards and design notes
- `cmake/` helper modules (e.g., BLAS, CUDA toggles)

## ParaView workflow (planned)

Examples will emit `.vtk` snapshots plus a `.pvd` collection file. Open `output.pvd` in ParaView and press Play to animate.

## Roadmap (phases)

1) Naive O(n²) N-body baseline
2) Optional cutoff prune
3) Spatial tree (quadtree/octree)
4) Dual-tree traversal scaffold
5) Exact cutoff via prune
6) Approximate long-range accept/refine (Barnes–Hut style)
7) Instrumentation and edge cases
8) CUDA translation using the same control structure

See `docs/CONVENTIONS.md` for coding style and design principles.
