#include "ramses/PoissonSolver.hpp"
#include "ramses/Parameters.hpp"
#include "ramses/Constants.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace ramses {

PoissonSolver::PoissonSolver(AmrGrid& grid) : grid_(grid) {
    phi.assign(grid_.ncell + 1, 0.0);
    rho.assign(grid_.ncell + 1, 0.0);
    f.allocate(grid_.ncell + 1, 3);
}

void PoissonSolver::solve(int ilevel) {
    if (grid_.count_grids_at_level(ilevel) == 0) return;

    if (params::verbose) {
        std::cout << "[PoissonSolver] Solving level " << ilevel << std::endl;
    }

    // 1. Prepare mask and BC-modified RHS
    make_fine_mask(ilevel);
    make_fine_bc_rhs(ilevel);

    // 2. Initial residual
    compute_residual(ilevel);
    real_t i_res_norm = compute_residual_norm(ilevel);

    // 3. Simple Iterative Solver (Gauss-Seidel)
    // For now, we implement a fixed number of iterations as a "simplified" solver
    const int MAXITER = 100;
    const real_t epsilon = 1e-4;
    
    real_t res_norm = i_res_norm;
    int iter = 0;
    while (iter < MAXITER && (iter == 0 || res_norm > epsilon * i_res_norm)) {
        smooth(ilevel, true);  // Red step
        smooth(ilevel, false); // Black step
        
        compute_residual(ilevel);
        res_norm = compute_residual_norm(ilevel);
        iter++;
        
        if (params::verbose && iter % 10 == 0) {
            std::cout << "  Iter " << iter << " Residual Norm: " << res_norm << std::endl;
        }
    }

    if (params::verbose) {
        std::cout << "  Finished at Iter " << iter << " Residual Norm: " << res_norm << std::endl;
    }
}

void PoissonSolver::smooth(int ilevel, bool redstep) {
    real_t dx = std::pow(0.5, ilevel);
    real_t dx2 = dx * dx;
    real_t dtwondim = static_cast<real_t>(constants::twondim);

    // Red/Black indices (1-based)
    std::vector<int> ired, iblack;
    if (grid_.ndim == 3) {
        ired = {1, 4, 6, 7};
        iblack = {2, 3, 5, 8};
    } else if (grid_.ndim == 2) {
        ired = {1, 4};
        iblack = {2, 3};
    } else {
        ired = {1};
        iblack = {2};
    }

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind0 = 0; ind0 < constants::twotondim / 2; ++ind0) {
                int ind = redstep ? ired[ind0] : iblack[ind0];
                int icell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;

                if (f(icell, 3) <= 0.0) continue;

                real_t nb_sum = 0.0;
                int igridn[7];
                grid_.get_nbor_grids(igrid, igridn);

                for (int idim = 0; idim < grid_.ndim; ++idim) {
                    for (int inbor = 0; inbor < 2; ++inbor) {
                        int ig = constants::iii[idim][inbor][ind - 1];
                        int ih = constants::jjj[idim][inbor][ind - 1];
                        int ig_nbor = igridn[ig];
                        if (ig_nbor > 0) {
                            int icell_nbor = grid_.ncoarse + (ih - 1) * grid_.ngridmax + ig_nbor;
                            nb_sum += phi[icell_nbor];
                        }
                    }
                }
                
                phi[icell] = (nb_sum - dx2 * f(icell, 2)) / dtwondim;
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::compute_residual(int ilevel) {
    real_t dx = std::pow(0.5, ilevel);
    real_t oneoverdx2 = 1.0 / (dx * dx);
    real_t dtwondim = static_cast<real_t>(constants::twondim);

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int icell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                if (f(icell, 3) <= 0.0) {
                    f(icell, 1) = 0.0;
                    continue;
                }

                real_t nb_sum = 0.0;
                int igridn[7];
                grid_.get_nbor_grids(igrid, igridn);

                for (int idim = 0; idim < grid_.ndim; ++idim) {
                    for (int inbor = 0; inbor < 2; ++inbor) {
                        int ig = constants::iii[idim][inbor][ind - 1];
                        int ih = constants::jjj[idim][inbor][ind - 1];
                        int ig_nbor = igridn[ig];
                        if (ig_nbor > 0) {
                            int icell_nbor = grid_.ncoarse + (ih - 1) * grid_.ngridmax + ig_nbor;
                            nb_sum += phi[icell_nbor];
                        }
                    }
                }
                
                // residual = (sum(phi_nb) - 2*ndim*phi) / dx^2 - RHS
                // We store -residual in f(icell, 1) to match Fortran
                f(icell, 1) = -oneoverdx2 * (nb_sum - dtwondim * phi[icell]) + f(icell, 2);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

real_t PoissonSolver::compute_residual_norm(int ilevel) {
    real_t norm2 = 0.0;
    int count = 0;
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int icell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                if (f(icell, 3) > 0.0) {
                    norm2 += f(icell, 1) * f(icell, 1);
                    count++;
                }
            }
            igrid = grid_.next[igrid - 1];
        }
    }
    return (count > 0) ? std::sqrt(norm2 / count) : 0.0;
}

void PoissonSolver::make_fine_mask(int ilevel) {
    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int icell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                f(icell, 3) = 1.0; // Active
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::make_fine_bc_rhs(int ilevel) {
    // Replicates make_fine_bc_rhs from Fortran
    real_t twopi = 2.0 * std::acos(-1.0);
    real_t scale = params::boxlen / static_cast<real_t>(params::nx); // Simplified nx_loc
    real_t fourpi = 2.0 * twopi * scale;
    
    // In cosmology mode: 1.5 * omega_m * aexp * scale
    // if (params::cosmo) fourpi = 1.5 * params::omega_m * params::aexp * scale;

    for (int icpu = 1; icpu <= grid_.ncpu; ++icpu) {
        int igrid = grid_.headl(icpu, ilevel);
        while (igrid > 0) {
            for (int ind = 1; ind <= constants::twotondim; ++ind) {
                int icell = grid_.ncoarse + (ind - 1) * grid_.ngridmax + igrid;
                f(icell, 2) = fourpi * (rho[icell] - grid_.rho_tot);
            }
            igrid = grid_.next[igrid - 1];
        }
    }
}

void PoissonSolver::vcycle(int ilevel) {
    // Placeholder for actual recursive V-cycle
    // Needs MG levels infrastructure
}

} // namespace ramses
