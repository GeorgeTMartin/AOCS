#include "simulation/initialization.hpp"
#include "sim_constants.hpp"
#include <json/reader.h>
#include <json/value.h>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace aocs::simulation {
namespace {

bool readVector3(
    const Json::Value& value,
    const std::string& name,
    Eigen::Vector3d& vector) {
    if (!value.isArray() || value.size() != 3 ||
        !value[0].isNumeric() || !value[1].isNumeric() || !value[2].isNumeric()) {
        std::cerr << name << " must be an array of three numbers\n";
        return false;
    }
    vector = Eigen::Vector3d(
        value[0].asDouble(),
        value[1].asDouble(),
        value[2].asDouble());
    if (!vector.allFinite()) {
        std::cerr << name << " must contain finite values\n";
        return false;
    }
    return true;
}

} // namespace

std::unique_ptr<Simulation> initializeSimulation() {
    std::ifstream configFile("config.json");
    if (!configFile) {
        configFile.clear();
        configFile.open("sim/config.json");
        if (!configFile) {
            std::cerr << "Unable to open config.json or sim/config.json\n";
            return nullptr;
        }
    }

    Json::Value config;
    Json::Reader reader;
    if (!reader.parse(configFile, config)) {
        std::cerr << reader.getFormattedErrorMessages();
        return nullptr;
    }

    auto simulation = std::make_unique<Simulation>();
    const Json::Value& initialConditions =
        config["simulation"]["initial_conditions"];
    if (!readVector3(
            initialConditions["r0_sc_eci_m"],
            "simulation.initial_conditions.r0_sc_eci_m",
            simulation->state.position_eci_m) ||
        !readVector3(
            initialConditions["v0_sc_eci_m"],
            "simulation.initial_conditions.v0_sc_eci_m",
            simulation->state.velocity_eci_m_s)) {
        return nullptr;
    }
    if (!initialConditions["t_J2000_s"].isNumeric()) {
        std::cerr << "simulation.initial_conditions.t_J2000_s must be seconds since J2000\n";
        return nullptr;
    }
    simulation->epoch_J2000_s = initialConditions["t_J2000_s"].asDouble();
    const Json::Value& simulationConfig = config["simulation"];
    if (!simulationConfig["orbit_count"].isUInt() ||
        simulationConfig["orbit_count"].asUInt() == 0 ||
        !simulationConfig["aocs_rate_hz"].isNumeric() ||
        !simulationConfig["gps_rate_hz"].isNumeric()) {
        std::cerr << "simulation.orbit_count, aocs_rate_hz, and gps_rate_hz must be positive values\n";
        return nullptr;
    }
    const double aocsRate_Hz = simulationConfig["aocs_rate_hz"].asDouble();
    const double gpsRate_Hz = simulationConfig["gps_rate_hz"].asDouble();
    simulation->mass_kg = config["spacecraft"]["mass_kg"].asDouble();
    if (!std::isfinite(simulation->epoch_J2000_s)) {
        std::cerr << "simulation.initial_conditions.t_J2000_s must be finite\n";
        return nullptr;
    }
    if (!std::isfinite(aocsRate_Hz) ||
        !std::isfinite(gpsRate_Hz) ||
        !std::isfinite(simulation->mass_kg) ||
        aocsRate_Hz <= 0.0 ||
        gpsRate_Hz <= 0.0 ||
        simulation->mass_kg <= 0.0) {
        std::cerr << "AOCS/GPS rates and spacecraft mass must be valid positive values\n";
        return nullptr;
    }
    const double gpsIntervalSteps = aocsRate_Hz / gpsRate_Hz;
    const double roundedGpsIntervalSteps = std::round(gpsIntervalSteps);
    if (roundedGpsIntervalSteps < 1.0 ||
        std::abs(gpsIntervalSteps - roundedGpsIntervalSteps) > 1e-9) {
        std::cerr << "AOCS rate must be an integer multiple of the GPS rate\n";
        return nullptr;
    }
    simulation->dt_sec = 1.0 / aocsRate_Hz;
    simulation->gpsSampleIntervalSteps =
        static_cast<std::uint64_t>(roundedGpsIntervalSteps);

    const double initialRadius_m = simulation->state.position_eci_m.norm();
    const double specificEnergy_m2_s2 =
        0.5 * simulation->state.velocity_eci_m_s.squaredNorm() -
        SimConstants::kEarthMu / initialRadius_m;
    if (!std::isfinite(initialRadius_m) || initialRadius_m <= 0.0 ||
        !std::isfinite(specificEnergy_m2_s2) || specificEnergy_m2_s2 >= 0.0) {
        std::cerr << "Initial conditions must define a bound Earth orbit\n";
        return nullptr;
    }
    const double semiMajorAxis_m =
        -SimConstants::kEarthMu / (2.0 * specificEnergy_m2_s2);
    const double orbitalPeriod_s =
        2.0 * SimConstants::kPi *
        std::sqrt(semiMajorAxis_m * semiMajorAxis_m * semiMajorAxis_m /
                  SimConstants::kEarthMu);
    const double requestedDuration_s =
        simulationConfig["orbit_count"].asDouble() * orbitalPeriod_s;
    simulation->simulationLength_s =
        std::round(requestedDuration_s / simulation->dt_sec) *
        simulation->dt_sec;

    const Json::Value& gravityConfig = config["environment"]["gravity"];
    if (!gravityConfig.isObject() || !gravityConfig["model"].isString()) {
        std::cerr << "environment.gravity.model must be 'point_mass' or 'egm2008'\n";
        return nullptr;
    }

    const std::string gravityModelName = gravityConfig["model"].asString();
    if (gravityModelName == "point_mass") {
        simulation->gravityModel = GravityModel::E_PointMass;
    } else if (gravityModelName == "egm2008") {
        simulation->gravityModel = GravityModel::E_EGM2008;
        if (!gravityConfig["degree"].isUInt() ||
            !gravityConfig["order"].isUInt()) {
            std::cerr << "environment.gravity.degree and order must be non-negative integers\n";
            return nullptr;
        }
        simulation->gravityDegree = gravityConfig["degree"].asUInt();
        simulation->gravityOrder = gravityConfig["order"].asUInt();
        try {
            simulation->egm2008 = std::make_unique<EGM2008Coefficients>(
                "data/earth_gravity/emg2008_coefficients.bin");
        } catch (const std::exception& error) {
            std::cerr << "Unable to load EGM2008 coefficients: "
                      << error.what() << '\n';
            return nullptr;
        }
        if (simulation->gravityDegree > simulation->egm2008->maxDegree ||
            simulation->gravityOrder > simulation->gravityDegree) {
            std::cerr << "EGM2008 degree must not exceed "
                      << simulation->egm2008->maxDegree
                      << " and order must not exceed degree\n";
            return nullptr;
        }
    } else {
        std::cerr << "Unsupported environment.gravity.model: "
                  << gravityModelName << '\n';
        return nullptr;
    }
    simulation->gravity =
        std::make_unique<Gravity>(simulation->egm2008.get());

    simulation->outputFile.open("output.csv");
    if (!simulation->outputFile) {
        std::cerr << "Unable to open output.csv for writing\n";
        return nullptr;
    }
    simulation->outputFile
        << "t_sec,r_x_m,r_y_m,r_z_m,v_x_m_s,v_y_m_s,v_z_m_s,"
        << "gps_time_s,gps_measurement_updated,"
        << "estimate_time_s,estimate_r_x_m,estimate_r_y_m,estimate_r_z_m,"
        << "estimate_v_x_m_s,estimate_v_y_m_s,estimate_v_z_m_s";
    for (std::size_t row = 0; row < 6; ++row) {
        for (std::size_t col = 0; col < 6; ++col) {
            simulation->outputFile << ",cov_" << row << '_' << col;
        }
    }
    simulation->outputFile << '\n' << std::setprecision(10);

    return simulation;
}

} // namespace aocs::simulation
