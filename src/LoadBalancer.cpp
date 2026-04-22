#include "ramses/LoadBalancer.hpp"
#include "ramses/MpiManager.hpp"
#include "ramses/Hilbert.hpp"
#include <iostream>
#include <algorithm>

namespace ramses {

void LoadBalancer::balance() {
    auto& mpi = MpiManager::instance();
    if (mpi.size() == 1) return;

    if (mpi.is_master()) {
        std::cout << "[LoadBalancer] Re-partitioning mesh..." << std::endl;
    }

    std::vector<int> cpu_map_new(grid_.ncell + 1);
    compute_new_cpu_map(cpu_map_new);
    
    move_grids(cpu_map_new);
}

void LoadBalancer::compute_new_cpu_map(std::vector<int>& cpu_map_new) {
    // 1. Generate Hilbert keys for all cells (simplified)
    // 2. Sort keys and split into N_CPU chunks with equal cell count (workload)
    
    // Placeholder: Split cells linearly
    auto& mpi = MpiManager::instance();
    int cells_per_cpu = grid_.ncell / mpi.size();
    
    for (int i = 1; i <= grid_.ncell; ++i) {
        int target_rank = (i - 1) / cells_per_cpu;
        cpu_map_new[i] = std::min(target_rank + 1, mpi.size());
    }
}

void LoadBalancer::move_grids(const std::vector<int>& cpu_map_new) {
    // This involves MPI_Isend/Irecv of oct data
    // Update local grid_.cpu_map
    for (int i = 1; i <= grid_.ncell; ++i) {
        grid_.cpu_map[i] = cpu_map_new[i];
    }
}

} // namespace ramses
