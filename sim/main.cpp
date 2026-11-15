#include "simulation/initialization.hpp"
#include "environment/earth_rotation/iau2000b.hpp"
#include "gnc_lib.hpp"
#include "ipc.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace {

constexpr double kPlaceholderGpsMeasurementVariance = 1.0;

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
    sensors.GPS_position_ECEF_m = {
        position_ecef_m.x(), position_ecef_m.y(), position_ecef_m.z()};
    const Eigen::Vector3d velocity_ecef_m_s =
        earthRotation.gcrsToItrs * simulation.state.velocity_eci_m_s -
        earthAngularVelocity_ecef_rad_s.cross(position_ecef_m);
    sensors.GPS_velocity_ECEF_m_s = {
        velocity_ecef_m_s.x(), velocity_ecef_m_s.y(), velocity_ecef_m_s.z()};
    sensors.GPS_time_s = t_J2000_s + aocs::ipc::kGpsSecondsAtJ2000;
    sensors.aocs_time_J2000_s = t_J2000_s;
    sensors.has_new_GPS_measurement = 1;
    for (std::size_t axis = 0; axis < aocs::ipc::kStateDimension; ++axis) {
        sensors.GPS_covariance_row_major[
            axis * aocs::ipc::kStateDimension + axis] =
            kPlaceholderGpsMeasurementVariance;
    }
    return sensors;
}

int runAOCSLoop(GncTransport& transport, const sensorMessage& sensors) {
    return transport.exchange(sensors);
}

void writeState(
    aocs::simulation::Simulation& simulation,
    const sensorMessage& sensors,
    const navigationEstimateMessage& estimate) {
    const Eigen::Vector3d& position = simulation.state.position_eci_m;
    const Eigen::Vector3d& velocity = simulation.state.velocity_eci_m_s;
    simulation.outputFile << simulation.elapsedTime_s << ','
                          << position.x() << ',' << position.y() << ',' << position.z()
                          << ',' << velocity.x() << ',' << velocity.y() << ','
                          << velocity.z() << ','
                          << sensors.GPS_time_s - aocs::ipc::kGpsSecondsAtJ2000 -
                                 simulation.epoch_J2000_s
                          << ',' << sensors.has_new_GPS_measurement
                          << ',' << estimate.time_J2000_s -
                                 simulation.epoch_J2000_s;
    for (double value : estimate.position_ECI_m) {
        simulation.outputFile << ',' << value;
    }
    for (double value : estimate.velocity_ECI_m_s) {
        simulation.outputFile << ',' << value;
    }
    for (double value : estimate.covariance_row_major) {
        simulation.outputFile << ',' << value;
    }
    simulation.outputFile << '\n';
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
        aocs::GncState gncState = aocs::initializeGnc();
        SharedMemoryGncTransport transport(
            sharedMemory.data(),
            [&gncState](sharedMemoryMap& memory) {
                return aocs::runGncStep(gncState, memory);
            });

        sensorMessage latestGpsMeasurement;
        std::uint64_t aocsStep = 0;
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
            if (aocsStep % simulation->gpsSampleIntervalSteps == 0) {
                latestGpsMeasurement = updateSensorModels(*simulation);
            } else {
                latestGpsMeasurement.has_new_GPS_measurement = 0;
                latestGpsMeasurement.aocs_time_J2000_s =
                    simulation->epoch_J2000_s + simulation->elapsedTime_s;
            }
            if (runAOCSLoop(transport, latestGpsMeasurement) != 0) {
                std::cerr << "GNC loop reported a failure\n";
                return 1;
            }
            const navigationEstimateMessage& navigationEstimate =
                sharedMemory.data().navigationEstimate;
            if (navigationEstimate.valid == 0) {
                std::cerr << "GNC did not publish a valid navigation estimate\n";
                return 1;
            }

            writeState(*simulation, latestGpsMeasurement, navigationEstimate);

            if (simulation->elapsedTime_s >= simulation->simulationLength_s) {
                std::cout << "Simulation ended at time: "
                          << simulation->elapsedTime_s << " sec\n";
                std::cout << "Final position: "
                          << simulation->state.position_eci_m.transpose() << "\n";
                std::cout << "Final velocity: "
                          << simulation->state.velocity_eci_m_s.transpose() << "\n";
                simulation->running = false;
            }
            ++aocsStep;
        }
    } catch (const std::exception& error) {
        std::cerr << "Simulation/GNC shared-memory failure: "
                  << error.what() << '\n';
        return 1;
    }

    return 0;
}
