#pragma once

#include "environment/earth_gravity/gravity.hpp"
#include "simulation/simulation.hpp"
#include <chrono>
#include <cstdint>
#include <fstream>
#include <memory>

namespace aocs::simulation {

struct Simulation {
    Propagation propagation;
    SpacecraftState state;
    std::chrono::steady_clock::time_point sim_start;
    double epoch_J2000_s = 0.0;
    double elapsedTime_s = 0.0;
    double simulationLength_s = 0.0;
    double dt_sec = 0.0;
    std::uint64_t gpsSampleIntervalSteps = 0;
    double mass_kg = 0.0;
    GravityModel gravityModel = GravityModel::E_PointMass;
    std::uint32_t gravityDegree = 0;
    std::uint32_t gravityOrder = 0;
    std::unique_ptr<EGM2008Coefficients> egm2008;
    std::unique_ptr<Gravity> gravity;
    std::ofstream outputFile;
    bool running = true;
};

std::unique_ptr<Simulation> initializeSimulation();

} // namespace aocs::simulation