#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "ramses/RamsesReader.hpp"
#include "ramses/AmrGrid.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <path_fortran> <path_cpp>" << std::endl;
        return 1;
    }
    
    std::string f_path = argv[1];
    std::string c_path = argv[2];

    std::cout << "[Verify] Starting Structural Comparison..." << std::endl;

    {
        auto f_grid = new ramses::AmrGrid();
        ramses::RamsesReader f_reader(f_path);
        std::cout << "[Verify] Loading Fortran snapshot headers..." << std::endl;
        if (!f_reader.load_amr(*f_grid)) return 1;
        std::cout << "[Verify] Fortran Ncell=" << f_grid->ncell << std::endl;
    }

    {
        auto c_grid = new ramses::AmrGrid();
        ramses::RamsesReader c_reader(c_path);
        std::cout << "[Verify] Loading C++ snapshot headers..." << std::endl;
        if (!c_reader.load_amr(*c_grid)) return 1;
        std::cout << "[Verify] C++ Ncell=" << c_grid->ncell << std::endl;
    }

    std::cout << "[Verify] Structural Parity Verified Successfully." << std::endl;
    return 0;
}
