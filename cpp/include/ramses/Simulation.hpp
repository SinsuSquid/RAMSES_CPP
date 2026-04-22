#ifndef RAMSES_SIMULATION_HPP
#define RAMSES_SIMULATION_HPP

#include "AmrGrid.hpp"
#include "HydroSolver.hpp"
#include "TreeUpdater.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Main simulation driver for RAMSES-CPP.
 * 
 * Orchestrates the AMR grid, Hydro solver, and Tree updates.
 */
class Simulation {
public:
    Simulation() : hydro_(grid_), updater_(grid_) {}

    /**
     * @brief Initializes the simulation from a namelist file.
     */
    void initialize(const std::string& nml_path);

    /**
     * @brief Executes the main time-stepping loop.
     */
    void run();

private:
    /**
     * @brief Performs a single AMR step at a specific level.
     * Replicates amr_step.f90 logic.
     */
    void amr_step(int ilevel);

    AmrGrid grid_;
    HydroSolver hydro_;
    TreeUpdater updater_;
    Config config_;

    real_t t_ = 0.0;
    real_t tend_ = 1.0;
    int nstep_ = 0;
    int nstepmax_ = 10;
};

} // namespace ramses

#endif // RAMSES_SIMULATION_HPP
