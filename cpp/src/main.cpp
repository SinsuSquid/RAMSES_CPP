#include <iostream>
#include "ramses/Simulation.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <namelist.nml>" << std::endl;
        return 1;
    }

    std::string nml_path = argv[1];
    
    ramses::Simulation sim;
    sim.initialize(nml_path);
    sim.run();

    return 0;
}
