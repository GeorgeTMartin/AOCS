#include "simulation/initialization.hpp"
#include "environment/earth_rotation/iau2000b.hpp"
#include "gnc_lib.hpp"
#include "ipc.hpp"
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
    const aocs::simulation::Simulation& simulation,
    const aocs::simulation::SpacecraftState& state,
    double elapsedTime_s) {
    const double t_J2000_s = simulation.epoch_J2000_s + elapsedTime_s;
    const aocs::environment::JulianDate date =
        aocs::environment::JulianDate::fromJ2000Seconds(t_J2000_s);
    const auto earthRotation = aocs::environment::EarthRotation(date, date).reduce();
    const Eigen::Vector3d position_ecef =
        earthRotation.gcrsToItrs * state.position_eci_m;
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

sensorMessage updateSensorModels(
    const aocs::simulation::Simulation& simulation) {
    const double t_J2000_s =
        simulation.epoch_J2000_s + simulation.elapsedTime_s;
    const aocs::environment::JulianDate date =
        aocs::environment::JulianDate::fromJ2000Seconds(t_J2000_s);
    const auto earthRotation =
        aocs::environment::EarthRotation(date, date).reduce();
    const Eigen::Vector3d position_ecef_m =
        earthRotation.gcrsToItrs * simulation.state.position_eci_m;
    const Eigen::Vector3d earthAngularVelocity_ecef_rad_s(
        0.0, 0.0, 7.2921150e-5);

    sensorMessage sensors;
    sensors.GPS_position_ECEF = position_ecef_m;
    sensors.GPS_velocity_ECEF =
        earthRotation.gcrsToItrs * simulation.state.velocity_eci_m_s -
        earthAngularVelocity_ecef_rad_s.cross(position_ecef_m);
    constexpr double kGpsSecondsAtJ2000 = 630763200.0;
    sensors.GPS_time_s = t_J2000_s + kGpsSecondsAtJ2000;
    return sensors;
}

int runAOCSLoop(GncTransport& transport, const sensorMessage& sensors) {
    return transport.exchange(sensors);
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

    try {
        SharedMemoryRegion sharedMemory;
        SharedMemoryGncTransport transport(
            sharedMemory.data(), &aocs::main);

        while (simulation->running) {
            const ActuatorState actuators = updateActuatorDynamics(*simulation);
            simulation->propagation.propagate(
                simulation->state,
                simulation->elapsedTime_s,
                simulation->dt_sec,
                [&](const aocs::simulation::SpacecraftState& state, double elapsedTime_s)
                    -> Eigen::Vector3d {
                    const EnvironmentalParameters environment =
                        fetchEnvironmentalParameters(*simulation, state, elapsedTime_s);
                    return Eigen::Vector3d(
                        (environment.gravityForce_eci_N + actuators.force_eci_N) /
                        simulation->mass_kg);
                });
            const sensorMessage sensors = updateSensorModels(*simulation);
            if (runAOCSLoop(transport, sensors) != 0) {
                std::cerr << "GNC loop reported a failure\n";
                return 1;
            }

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
    } catch (const std::exception& error) {
        std::cerr << "Simulation/GNC shared-memory failure: "
                  << error.what() << '\n';
        return 1;
    }

    return 0;
}
