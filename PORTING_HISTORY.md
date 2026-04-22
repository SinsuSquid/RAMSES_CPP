# RAMSES-2025 C++ Porting History

This file tracks the architectural decisions and progress of the port from Fortran to C++.

## Phase 0: Infrastructure & Foundation

### [2026-04-22] - Project Initialization
- **C++ Environment:** Created `cpp/` directory structure and `CMakeLists.txt`.
- **Standards:** Targeting C++17.
- **Type Mapping:** Created `Types.hpp` for bit-perfect type alignment.
- **Constants:** Ported foundation into `Constants.hpp` and `Parameters.hpp`.

## Phase 4: AMR Tree Logic & Refinement (Current)

### [2026-04-22] - Tree Traversal & Coarse Refinement
- **Neighbor Search:** Implemented `get_nbor_grids` and `get_nbor_cells` in `AmrGrid` using dimension-agnostic lookup tables.
- **Refinement Logic:** Implemented `TreeUpdater::refine_coarse` with periodic boundary support.
- **Verification:** Successfully verified coarse-level oct creation and neighbor connectivity.

## Phase 5: Validation Infrastructure (Current)

### [2026-04-22] - Fortran Binary Bridge
- **RamsesReader:** Implemented a C++ utility to read RAMSES unformatted Fortran binary files.
- **Snapshot Loading:** Developed `load_amr` and `load_hydro` to reconstruct the full state (Grid + Conservative Physics) from Fortran output.
- **Physics Conversion:** Implemented automatic conversion from primitive variables (saved in RAMSES) to conservative variables (used in C++ `uold`).
- **Parity:** This enables bit-for-bit comparison between the legacy Fortran code and the new C++ port.

## Phase 6: Hydro Solver (Current)

### [2026-04-22] - Godunov Solver & Riemann Physics
- **HydroSolver Class:** Main wrapper for hydrodynamics and state management.
- **RiemannSolver:** Standalone LLF solver implemented and verified.
- **SlopeLimiter:** Implemented bit-perfect TVD slope limiters (MinMod, MonCen, van Leer).
- **MUSCL Tracing:** Implemented MUSCL-Hancock prediction logic in `Muscl` class for interface state reconstruction.
- **Verification:** Verified flux calculations and slope limiting against standard test cases.

### Architectural Decisions
1. **Physics Modularity:** Slope limiting, Riemann solving, and MUSCL tracing are implemented as independent, stateless classes to facilitate rigorous unit testing.
2. **Stencil Handling:** Initial framework for `godfine1` batch processing of grids is established, matching RAMSES' vector-sweep optimization.

## Upcoming Milestones
- [ ] **MUSCL Reconstruction:** Implement slope limiters (MinMod, Moncen).
- [ ] **Unsplit Scheme:** Integrate the 3D unsplit Godunov solver.
- [ ] **Source Terms:** Add gravity and PdV source terms.
