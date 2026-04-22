#include <iostream>
#include "ramses/Simulation.hpp"
#include "ramses/MpiManager.hpp"

int main(int argc, char** argv) {
    auto& mpi = ramses::MpiManager::instance();
    mpi.init(argc, argv);

    if (argc < 2) {
        if (mpi.is_master()) {
            std::cout << "Usage: " << argv[0] << " <namelist.nml>" << std::endl;
        }
        mpi.finalize();
        return 1;
    }

    std::string nml_path = argv[1];
    
    ramses::Simulation sim;
    sim.initialize(nml_path);
    sim.run();

    mpi.finalize();
    return 0;
}
