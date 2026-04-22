#include "ramses/ParticleSolver.hpp"
#include <iostream>
#include <cmath>

namespace ramses {

void ParticleSolver::move_fine(int ilevel, real_t dt) {
    // In production, we'd loop over active grids at ilevel.
    // For now, we traverse the entire grid system to find particles.
    std::vector<int> ind_part;
    
    // Simplification: Loop over all particles and move those belonging to ilevel
    for (int ipart = 1; ipart <= ps_.npartmax; ++ipart) {
        if (ps_.idp[ipart - 1] > 0 && ps_.levelp[ipart - 1] == ilevel) {
            ind_part.push_back(ipart);
            if (ind_part.size() >= 100) { // Batch processing
                move_particles(ind_part, dt);
                ind_part.clear();
            }
        }
    }
    if (!ind_part.empty()) move_particles(ind_part, dt);
}

void ParticleSolver::move_particles(const std::vector<int>& ind_part, real_t dt) {
    const real_t boxlen = 1.0; // Should come from params

    for (int ipart : ind_part) {
        // Simple advection: x = x + v * dt
        for (int idim = 1; idim <= NDIM; ++idim) {
            real_t& x = ps_.get_xp(ipart, idim);
            real_t v = ps_.get_vp(ipart, idim);
            x += v * dt;

            // Periodic Boundaries
            if (x < 0.0) x += boxlen;
            if (x >= boxlen) x -= boxlen;
        }
    }
}

} // namespace ramses
