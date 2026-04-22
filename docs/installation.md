---
layout: default
---

# Installation

## Prerequisites
- **CMake** (>= 3.15)
- **C++17 compliant compiler** (e.g., GCC 9+, Clang)
- **MPI** (Optional)

## Building from Source

1. Clone the repository:
   ```bash
   git clone https://github.com/SinsuSquid/RAMSES_CPP.git
   cd RAMSES_CPP
   ```

2. Build:
   ```bash
   mkdir build && cd build
   cmake ..
   make -j$(nproc)
   ```
