#ifndef RAMSES_LOAD_BALANCER_HPP
#define RAMSES_LOAD_BALANCER_HPP

#include "AmrGrid.hpp"
#include "ParticleSystem.hpp"
#include <vector>

namespace ramses {

/**
 * @brief Handles MPI domain decomposition and load balancing.
 * 
 * Replicates logic from amr/load_balance.f90.
 */
class LoadBalancer {
public:
    LoadBalancer(AmrGrid& grid, ParticleSystem& ps) : grid_(grid), ps_(ps) {}

    /**
     * @brief Performs full load balancing across all MPI ranks.
     */
    void balance();

private:
    /**
     * @brief Computes a new CPU map based on Hilbert curve ordering.
     */
    void compute_new_cpu_map(std::vector<int>& cpu_map_new);

    /**
     * @brief Physically moves grids between MPI ranks.
     */
    void move_grids(const std::vector<int>& cpu_map_new);

    AmrGrid& grid_;
    ParticleSystem& ps_;
};

} // namespace ramses

#endif // RAMSES_LOAD_BALANCER_HPP
