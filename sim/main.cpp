#include "simulation/initialization.hpp"
#include "environment/earth_rotation/iau2000b.hpp"
#include "astro/ijk_keplerian.hpp"
#include "gnc_lib.hpp"
#include "ipc.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

namespace {

constexpr double kPlaceholderGpsMeasurementVariance = 1.0;
constexpr std::size_t kTelemetryLineCount = 7;

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

std::string renderTelemetry(
    const aocs::simulation::Simulation& simulation,
    bool useColor) {
    const orbitalElements elements = IJKKeplerian::convert_ijk_to_keplerian(
        simulation.state.position_eci_m,
        simulation.state.velocity_eci_m_s);
    const double progress = std::clamp(
        simulation.elapsedTime_s / simulation.simulationLength_s,
        0.0,
        1.0);
    constexpr int progressBarWidth = 30;
    const int progressBarFill =
        static_cast<int>(std::round(progress * progressBarWidth));
    constexpr double kRadiansToDegrees = 180.0 / 3.14159265358979323846;
    const char* blue = useColor ? "\033[38;2;74;117;158m" : "";
    const char* orange = useColor ? "\033[38;2;204;132;82m" : "";
    const char* reset = useColor ? "\033[0m" : "";

    std::ostringstream display;
    display << std::fixed << std::setprecision(3)
            << blue << "SIMULATION TELEMETRY" << reset << '\n'
            << blue << "Time: " << reset << orange << simulation.elapsedTime_s
            << reset << blue << " / " << reset << orange
            << simulation.simulationLength_s << reset << blue
            << " s | J2000: " << reset << orange
            << simulation.epoch_J2000_s + simulation.elapsedTime_s << reset
            << blue << " s" << reset << '\n'
            << std::setprecision(2)
            << blue << "Orbit: a = " << reset << orange << elements.sma << reset
            << blue << " m  e = " << reset << orange << elements.e << reset
            << blue << "  i = " << reset << orange
            << elements.i * kRadiansToDegrees << reset << blue << " deg"
            << reset << '\n'
            << blue << "       RAAN = " << reset << orange
            << elements.omega * kRadiansToDegrees << reset << blue
            << " deg  argp = " << reset << orange
            << elements.argp * kRadiansToDegrees << reset << blue
            << " deg  nu = " << reset << orange
            << elements.nu * kRadiansToDegrees << reset << blue << " deg"
            << reset << '\n'
            << std::setprecision(1)
            << blue << std::left << std::setw(20) << "ECI position [m]:"
            << reset << orange << std::right
            << std::setw(14) << simulation.state.position_eci_m.x()
            << std::setw(14) << simulation.state.position_eci_m.y()
            << std::setw(14) << simulation.state.position_eci_m.z()
            << reset << '\n'
            << blue << std::left << std::setw(20) << "ECI velocity [m/s]:"
            << reset << orange << std::right
            << std::setw(14) << simulation.state.velocity_eci_m_s.x()
            << std::setw(14) << simulation.state.velocity_eci_m_s.y()
            << std::setw(14) << simulation.state.velocity_eci_m_s.z()
            << reset << '\n'
            << blue << '[' << reset << orange
            << std::string(progressBarFill, '=') << reset << blue
            << std::string(progressBarWidth - progressBarFill, ' ') << ']'
            << reset << ' ' << orange << std::setprecision(1) << progress * 100.0
            << reset << blue << "%  Simulation Rate: " << reset << orange
            << simulation.elapsedTime_s /
                   std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - simulation.sim_start)
                       .count()
            << reset << blue << "x Real-Time" << reset;
    return display.str();
}

void publishTelemetry(
    const aocs::simulation::Simulation& simulation,
    bool interactive,
    bool& hasDrawnTelemetry) {
    if (!interactive && simulation.running) {
        return;
    }

    if (interactive && hasDrawnTelemetry) {
        std::cout << "\033[" << kTelemetryLineCount << 'A';
    }
    std::istringstream lines(renderTelemetry(simulation, interactive));
    std::string line;
    while (std::getline(lines, line)) {
        if (interactive) {
            std::cout << "\r\033[2K";
        }
        std::cout << line << '\n';
    }
    std::cout.flush();
    hasDrawnTelemetry = true;
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
        aocs::GNCExecutive gncExecutive = aocs::initializeGnc();
        SharedMemoryGncTransport transport(
            sharedMemory.data(),
            [&gncExecutive](sharedMemoryMap& memory) {
                return aocs::runGncStep(gncExecutive, memory);
            });

        const bool interactive = isatty(STDOUT_FILENO) != 0;
        bool hasDrawnTelemetry = false;
        simulation->sim_start = std::chrono::steady_clock::now();
        auto lastTelemetryUpdate = simulation->sim_start;
        publishTelemetry(*simulation, interactive, hasDrawnTelemetry);
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

            const auto now = std::chrono::steady_clock::now();
            if (now - lastTelemetryUpdate >= std::chrono::milliseconds(100)) {
                publishTelemetry(*simulation, interactive, hasDrawnTelemetry);
                lastTelemetryUpdate = now;
            }
            if (simulation->elapsedTime_s >= simulation->simulationLength_s) {
                simulation->running = false;
            }
            ++aocsStep;
        }
        publishTelemetry(*simulation, interactive, hasDrawnTelemetry);
    } catch (const std::exception& error) {
        std::cerr << "Simulation/GNC shared-memory failure: "
                  << error.what() << '\n';
        return 1;
    }

    return 0;
}
