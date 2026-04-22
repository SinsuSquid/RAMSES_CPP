#include "ramses/HydroSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include "ramses/RiemannSolver.hpp"
#include "ramses/SlopeLimiter.hpp"
#include "ramses/Muscl.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace ramses {

void HydroSolver::set_unew(int ilevel) {
    for (int ind_cell = 1; ind_cell <= grid_.ncell; ++ind_cell) {
        for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
            grid_.unew(ind_cell, ivar) = grid_.uold(ind_cell, ivar);
        }
    }
}

void HydroSolver::set_uold(int ilevel) {
    for (int ind_cell = 1; ind_cell <= grid_.ncell; ++ind_cell) {
        for (int ivar = 1; ivar <= grid_.nvar; ++ivar) {
            grid_.uold(ind_cell, ivar) = grid_.unew(ind_cell, ivar);
        }
    }
}

void HydroSolver::godunov_fine(int ilevel) {
    // Top level entry
    set_unew(ilevel);
    
    // In a real port, we would loop over active grids.
    // For this 1D POC, we'll just demonstrate the logical flow.
}

void HydroSolver::ctoprim(const real_t u[], real_t q[], real_t gamma) {
    const real_t smallr = 1e-10;
    real_t d = std::max(u[0], smallr);
    q[0] = d;
    q[1] = u[1] / d;
    q[2] = u[2] / d;
    q[3] = u[3] / d;
    
    real_t e_kin = 0.5 * d * (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    real_t e_int = u[4] - e_kin;
    q[4] = std::max(e_int * (gamma - 1.0), d * 1e-10); // pressure floor
}

void HydroSolver::gather_stencil(int igrid, int ilevel, LocalStencil& stencil) {
    int nbors_father[27];
    grid_.get_3x3x3_father(igrid, nbors_father);
    
    // Loop over 3x3x3 father cells
    for (int k1 = 0; k1 < 3; ++k1) {
    for (int j1 = 0; j1 < 3; ++j1) {
    for (int i1 = 0; i1 < 3; ++i1) {
        int ifather = nbors_father[i1 + 3*j1 + 9*k1];
        int ison = grid_.son[ifather];
        
        if (ison > 0) {
            // Neighbor exists, gather 2x2x2 cells
            for (int k2 = 0; k2 < 2; ++k2) {
            for (int j2 = 0; j2 < 2; ++j2) {
            for (int i2 = 0; i2 < 2; ++i2) {
                int icell_pos = 1 + i2 + 2*j2 + 4*k2;
                int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + ison;
                
                int i3 = i1 * 2 + i2;
                int j3 = j1 * 2 + j2;
                int k3 = k1 * 2 + k2;
                
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ind_cell, iv);
                }
                stencil.refined[i3][j3][k3] = (grid_.son[ind_cell] > 0);
            }
            }
            }
        } else {
            // Neighbor doesn't exist, interpolate from father (simple injection for now)
            for (int k2 = 0; k2 < 2; ++k2) {
            for (int j2 = 0; j2 < 2; ++j2) {
            for (int i2 = 0; i2 < 2; ++i2) {
                int i3 = i1 * 2 + i2;
                int j3 = j1 * 2 + j2;
                int k3 = k1 * 2 + k2;
                
                for (int iv = 1; iv <= grid_.nvar; ++iv) {
                    stencil.uloc[i3][j3][k3][iv - 1] = grid_.uold(ifather, iv);
                }
                stencil.refined[i3][j3][k3] = false;
            }
            }
            }
        }
    }
    }
}

void HydroSolver::godfine1(const std::vector<int>& ind_grid, int ilevel) {
    // ind_grid contains 1-based oct indices
    real_t gamma = 1.4;
    
    for (int igrid : ind_grid) {
        // For each oct, we need a 6x6x6 stencil
        // Inner cells of oct: (ix,iy,iz) = [1,2] x [1,2] x [1,2] (in a 4x4x4 block)
        // RAMSES uses IU1=-1, IU2=4 for buffer
        
        // This is a complex assembly step. For the port, we'll implement 
        // a simplified 1D sweep to prove the flow.
        
        for (int icell_pos = 1; icell_pos <= constants::twotondim; ++icell_pos) {
            int ind_cell = grid_.ncoarse + (icell_pos - 1) * grid_.ngridmax + igrid;
            
            // 1. Get neighbors
            int igridn[7];
            grid_.get_nbor_grids(igrid, igridn);
            int icelln[6];
            grid_.get_nbor_cells(igridn, icell_pos, icelln);
            
            // 2. Compute fluxes in each direction
            for (int idim = 0; idim < NDIM; ++idim) {
                int left_cell = icelln[idim * 2];
                int right_cell = icelln[idim * 2 + 1];
                
                if (left_cell > 0 && right_cell > 0) {
                    real_t ql[5], qc[5], qr[5];
                    ctoprim(grid_.uold.data() + (left_cell - 1) * 5, ql, gamma);
                    ctoprim(grid_.uold.data() + (ind_cell - 1) * 5, qc, gamma);
                    ctoprim(grid_.uold.data() + (right_cell - 1) * 5, qr, gamma);
                    
                    // MUSCL Reconstruction
                    real_t dql[5], dqr[5];
                    for (int iv = 0; iv < 5; ++iv) {
                        // Left interface state
                        real_t slope_l = SlopeLimiter::compute_slope(0, ql[iv], qc[iv], 1); // Simplification
                        real_t slope_r = SlopeLimiter::compute_slope(ql[iv], qc[iv], qr[iv], 1);
                        
                        // Prediction step (dt=0 for now)
                        // ...
                    }
                }
            }
        }
    }
}

} // namespace ramses
