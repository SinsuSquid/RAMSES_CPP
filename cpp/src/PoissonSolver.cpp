#include "ramses/PoissonSolver.hpp"
#include <iostream>

namespace ramses {

void PoissonSolver::solve(int ilevel) {
    if (grid_.ncoarse == 0) return;
    
    // Placeholder for Multigrid cycle:
    smooth(ilevel);
}

void PoissonSolver::smooth(int ilevel) {
    // Gauss-Seidel placeholder
}

} // namespace ramses
