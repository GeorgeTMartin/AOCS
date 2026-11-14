#include "simulation/initialization.hpp"
#include "environment/earth_rotation/iau2000b.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace {

struct EnvironmentalParameters {
    Eigen::Vector3d gravityForce_eci_N;
};

struct ActuatorState {
    Eigen::Vector3d force_eci_N = Eigen::Vector3d::Zero();
};

EnvironmentalParameters fetchEnvironmentalParameters(
    const aocs::simulation::Simulation& simulation) {
    const double t_J2000_s = simulation.epoch_J2000_s + simulation.elapsedTime_s;
    const aocs::environment::JulianDate date =
        aocs::environment::JulianDate::fromJ2000Seconds(t_J2000_s);
    const auto earthRotation = aocs::environment::EarthRotation(date, date).reduce();
    const Eigen::Vector3d position_ecef =
        earthRotation.gcrsToItrs * simulation.state.position_eci_m;
    const Eigen::Vector3d gravityForce_ecef = simulation.gravity->computeGravity(
        position_ecef,
        simulation.mass_kg,
        simulation.gravityModel,
        simulation.gravityDegree,
        simulation.gravityOrder);

    return {earthRotation.gcrsToItrs.transpose() * gravityForce_ecef};
}

ActuatorState updateActuatorDynamics(const aocs::simulation::Simulation&) {
    return {};
}

void propagateState(
    aocs::simulation::Simulation& simulation,
    const EnvironmentalParameters& environment,
    const ActuatorState& actuators) {
    const Eigen::Vector3d acceleration_eci_m_s2 =
        (environment.gravityForce_eci_N + actuators.force_eci_N) / simulation.mass_kg;
    simulation.state.velocity_eci_m_s += acceleration_eci_m_s2 * simulation.dt_sec;
    simulation.state.position_eci_m +=
        simulation.state.velocity_eci_m_s * simulation.dt_sec;
    simulation.elapsedTime_s += simulation.dt_sec;
}

void updateSensorModels(const aocs::simulation::Simulation&) {
    // Sensor models are not configured yet.
}

void runAOCSLoop(aocs::simulation::Simulation&) {
    // The AOCS interface is not connected yet.
}

void writeState(aocs::simulation::Simulation& simulation) {
    const Eigen::Vector3d& position = simulation.state.position_eci_m;
    const Eigen::Vector3d& velocity = simulation.state.velocity_eci_m_s;
    simulation.outputFile << simulation.elapsedTime_s << ','
                          << position.x() << ',' << position.y() << ',' << position.z()
                          << ',' << velocity.x() << ',' << velocity.y() << ',' << velocity.z()
                          << '\n';
}

} // namespace

int main() {
    std::unique_ptr<aocs::simulation::Simulation> simulation =
        aocs::simulation::initializeSimulation();
    if (!simulation) {
        return 1;
    }

    while (simulation->running) {
        const EnvironmentalParameters environment =
            fetchEnvironmentalParameters(*simulation);
        const ActuatorState actuators = updateActuatorDynamics(*simulation);
        propagateState(*simulation, environment, actuators);
        updateSensorModels(*simulation);
        runAOCSLoop(*simulation);

        writeState(*simulation);

        if (simulation->elapsedTime_s >= simulation->simulationLength_s) {
            std::cout << "Simulation ended at time: "
                      << simulation->elapsedTime_s << " sec\n";
            std::cout << "Final position: "
                      << simulation->state.position_eci_m.transpose() << "\n";
            std::cout << "Final velocity: "
                      << simulation->state.velocity_eci_m_s.transpose() << "\n";
            simulation->running = false;
        }
    }

    return 0;
}
