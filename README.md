# AOCS

This repository is intended to provide the basis of development for primitive AOCS capabilities, with an associated simulation framework.

## Build

Requirements: CMake 3.20 or newer, a C++17 compiler, and the JsonCpp development package. Eigen is included in the repository.

From a fresh clone, configure and build with:

```sh
cmake -S . -B build
cmake --build build
```

The generated build directory is local and is not tracked by Git.
