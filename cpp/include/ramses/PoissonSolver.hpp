#ifndef RAMSES_POISSON_SOLVER_HPP
#define RAMSES_POISSON_SOLVER_HPP

#include "AmrGrid.hpp"

namespace ramses {

/**
 * @brief Multigrid-based Poisson solver for gravity.
 * 
 * Replicates logic from poisson/ directory.
 */
class PoissonSolver {
public:
    PoissonSolver(AmrGrid& grid) : grid_(grid) {}

    /**
     * @brief Solves Del^2 Phi = 4*pi*G*rho.
     */
    void solve(int ilevel);

private:
    AmrGrid& grid_;
    
    /**
     * @brief Performs a Gauss-Seidel smoothing step.
     */
    void smooth(int ilevel);
};

} // namespace ramses

#endif // RAMSES_POISSON_SOLVER_HPP
