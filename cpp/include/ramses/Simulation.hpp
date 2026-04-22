#ifndef RAMSES_SIMULATION_HPP
#define RAMSES_SIMULATION_HPP

#include "AmrGrid.hpp"
#include "HydroSolver.hpp"
#include "PoissonSolver.hpp"
#include "TreeUpdater.hpp"
#include "Config.hpp"

namespace ramses {

/**
 * @brief Main simulation driver for RAMSES-CPP.
 */
class Simulation {
public:
    Simulation() : hydro_(grid_), poisson_(grid_), updater_(grid_) {}

    void initialize(const std::string& nml_path);
    void run();

private:
    void amr_step(int ilevel);
    void dump_snapshot(int iout);

    AmrGrid grid_;
    HydroSolver hydro_;
    PoissonSolver poisson_;
    TreeUpdater updater_;
    Config config_;

    real_t t_ = 0.0;
    real_t tend_ = 1.0;
    int nstep_ = 0;
    int nstepmax_ = 10;
};

} // namespace ramses

#endif // RAMSES_SIMULATION_HPP
