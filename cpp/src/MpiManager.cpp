#include "ramses/MpiManager.hpp"
#include <iostream>

namespace ramses {

void MpiManager::init(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    
    if (is_master()) {
        std::cout << "[MPI] Initialized with " << size_ << " ranks." << std::endl;
    }
}

void MpiManager::finalize() {
    MPI_Finalize();
}

} // namespace ramses
