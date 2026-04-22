#ifndef RAMSES_POISSON_SOLVER_HPP
#define RAMSES_POISSON_SOLVER_HPP

#include "AmrGrid.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Multigrid-based Poisson solver for gravity.
 * 
 * Replicates logic from poisson/ directory.
 */
class PoissonSolver {
public:
    PoissonSolver(AmrGrid& grid);

    /**
     * @brief Solves Del^2 Phi = 4*pi*G*rho.
     */
    void solve(int ilevel);

    // Fields
    std::vector<real_t> phi;      // Potential
    std::vector<real_t> rho;      // Density (RHS)
    Field<real_t> f;              // 1: residual, 2: BC-modified RHS, 3: mask

private:
    AmrGrid& grid_;
    
    /**
     * @brief Performs a Gauss-Seidel smoothing step.
     */
    void smooth(int ilevel, bool redstep);

    /**
     * @brief Computes the residual Del^2 Phi - 4*pi*G*rho.
     */
    void compute_residual(int ilevel);

    /**
     * @brief Performs a V-cycle step (recursive).
     */
    void vcycle(int ilevel);

    // Helper methods
    void make_fine_mask(int ilevel);
    void make_fine_bc_rhs(int ilevel);
    real_t compute_residual_norm(int ilevel);
};

} // namespace ramses

#endif // RAMSES_POISSON_SOLVER_HPP
