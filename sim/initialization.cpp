#include "simulation/initialization.hpp"
#include <json/reader.h>
#include <json/value.h>
#include <cmath>
#include <exception>
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
    simulation->simulationLength_s = config["simulation"]["sim_length_sec"].asDouble();
    simulation->dt_sec = config["simulation"]["dt_sec"].asDouble();
    simulation->mass_kg = config["spacecraft"]["mass_kg"].asDouble();
    if (!std::isfinite(simulation->epoch_J2000_s)) {
        std::cerr << "simulation.initial_conditions.t_J2000_s must be finite\n";
        return nullptr;
    }
    if (!std::isfinite(simulation->simulationLength_s) ||
        !std::isfinite(simulation->dt_sec) ||
        !std::isfinite(simulation->mass_kg) ||
        simulation->simulationLength_s <= 0.0 ||
        simulation->dt_sec <= 0.0 ||
        simulation->mass_kg <= 0.0) {
        std::cerr << "Simulation length, timestep, and spacecraft mass must be valid positive values\n";
        return nullptr;
    }

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
        << "t_sec,r_x_m,r_y_m,r_z_m,v_x_m_s,v_y_m_s,v_z_m_s\n";

    return simulation;
}

} // namespace aocs::simulation
