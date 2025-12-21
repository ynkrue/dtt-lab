# Coding conventions

- Namespaces: `dtt::…`; files `snake_case`; types `PascalCase`; functions/vars `snake_case`; constants `kCamelCase`.
- Headers: `#pragma once`; declarations in headers, definitions in `src`; include order = own header, standard library, third-party.
- Interfaces: prefer `std::span`/`std::string_view`; pass by `const&`; mark everything `const` you can.
- Memory: RAII; no raw `new`/`delete`; avoid owning globals; deterministic destruction order.
- Errors: avoid exceptions for control flow; return status/enum or optional/expected style; log via `std::cerr` when needed.
- Layouts: favor SoA for particle data and explicit strides; avoid virtual dispatch in hot paths; keep structs trivially copyable when possible.
- Parallel/GPU prep: avoid recursion in traversal APIs; expose explicit worklists; keep data structures host/device friendly where reasonable.
- Testing: GoogleTest; deterministic seeds; keep benchmarks in `benchmarks/` not in tests.
- Formatting: `.clang-format` (LLVM base, 4 spaces, 100 cols, sorted includes); `clang-tidy` optional; treat warnings as errors in library/tests.
- VTK output: emit legacy `.vtk` or `.pvd` collections from examples; keep I/O helpers small and dependency-free.
- BLAS: link via `find_package(BLAS REQUIRED)`; wrap calls behind thin helpers so we can swap providers (OpenBLAS, Accelerate, MKL).
