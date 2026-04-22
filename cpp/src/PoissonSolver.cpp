#include "ramses/PoissonSolver.hpp"
#include <iostream>

namespace ramses {

void PoissonSolver::solve(int ilevel) {
    if (grid_.ncoarse == 0) return;
    
    std::cout << "[Poisson] Solving gravity for level " << ilevel << "..." << std::endl;
    
    // Placeholder for Multigrid cycle:
    // 1. Smooth (Gauss-Seidel)
    // 2. Compute residual
    // 3. Restrict to coarser level
    // 4. Solve coarse problem
    // 5. Prolongate error and correct
    // 6. Smooth again
    
    smooth(ilevel);
}

void PoissonSolver::smooth(int ilevel) {
    // Gauss-Seidel placeholder
}

} // namespace ramses
