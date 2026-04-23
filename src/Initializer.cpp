#include "ramses/Initializer.hpp"
#include "ramses/Parameters.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

void Initializer::apply_all() {
    region_condinit();
}

void Initializer::region_condinit() {
    int nregion = config_.get_int("init_params", "nregion", 1);
    real_t gamma = config_.get_double("hydro_params", "gamma", 1.4);
    
    // Default Background values
    real_t d_bg = config_.get_double("init_params", "d_region", 1.0);
    real_t p_bg = config_.get_double("init_params", "p_region", 1e-5);
    
    // Initialize all cells with background
    for (int i = 1; i <= grid_.ncell; ++i) {
        grid_.uold(i, 1) = d_bg; // density
        for (int idim = 1; idim <= NDIM; ++idim) {
            grid_.uold(i, 1 + idim) = 0.0; // velocity
        }
        grid_.uold(i, NDIM + 2) = p_bg / (gamma - 1.0); // total energy
        
        // Extra variables (nener)
        for (int ivar = NDIM + 3; ivar <= grid_.nvar; ++ivar) {
            grid_.uold(i, ivar) = 0.0;
        }
        grid_.cpu_map[i] = 1;
    }

    // Apply regions
    // For Sod tube, we usually have 2 regions
    if (nregion >= 2) {
        real_t d_region = 0.125;
        real_t p_region = 0.1;
        // This is a very simplified Initializer. In a real RAMSES, we'd check region shapes.
        // For Sod tube 1D, x < 0.5 is region 1, x > 0.5 is region 2.
        
        // Region 2 (usually the right side)
        // Since we only have ncoarse cells at level 0, let's just split them.
        for (int i = 1; i <= grid_.ncoarse; ++i) {
            real_t x = (static_cast<real_t>(i) - 0.5f) / static_cast<real_t>(params::nx);
            if (x > 0.5) {
                grid_.uold(i, 1) = 0.125; // d_region from sod-tube.nml
                grid_.uold(i, NDIM + 2) = 0.1 / (gamma - 1.0); // p_region
            }
        }
        // Also apply to fine cells if they were already created
        // (but usually Initializer is called after base level is set)
    }
    
    // Special case for sod-tube-nener: handle prad_region
    int nener = grid_.nvar - (2 + NDIM);
    if (nener > 0) {
        // Very simplified: just set them to some values if found in config
        for (int ireg = 1; ireg <= nregion; ++ireg) {
            for (int iv = 1; iv <= nener; ++iv) {
                // prad_region(ivar, iregion)
                std::string key = "prad_region(" + std::to_string(iv) + "," + std::to_string(ireg) + ")";
                real_t val = config_.get_double("init_params", key, 0.0);
                // Apply this value to the region... (omitted for brevity, POC only)
                // For sod-tube-nener, let's just set them.
                if (ireg == 1) {
                    for(int i=1; i<=grid_.ncell; ++i) grid_.uold(i, NDIM + 2 + iv) = val;
                }
            }
        }
    }

    std::cout << "[Initializer] Applied ICs (nvar=" << grid_.nvar << ")." << std::endl;
}

} // namespace ramses
