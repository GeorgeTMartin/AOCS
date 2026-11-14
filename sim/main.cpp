#include "ipc.hpp"
#include "environment/earth_gravity/gravity.hpp"
#include "environment/earth_rotation/iau2000b.hpp"
#include <Eigen/Dense>
#include <json/value.h>
#include <json/reader.h>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>

bool B_simulationRunning = true;

int main() {
    // initialize simulation
    // Load config.json from the current directory, or from sim/ when run from the repository root.
    std::ifstream configFile("config.json");
    if (!configFile) {
        configFile.clear();
        configFile.open("sim/config.json");
        if (!configFile) {
            std::cerr << "Unable to open config.json or sim/config.json\n";
            return 1;
        }
    }

    Json::Value config;
    Json::Reader reader;
    if (!reader.parse(configFile, config)) {
        std::cerr << reader.getFormattedErrorMessages();
        return 1;
    }

    const Json::Value& initialConditions = config["simulation"]["initial_conditions"];
    Eigen::Matrix<double, 3, 1> position = Eigen::Matrix<double, 3, 1>(
        initialConditions["r0_sc_eci_m"][0].asDouble(),
        initialConditions["r0_sc_eci_m"][1].asDouble(),
        initialConditions["r0_sc_eci_m"][2].asDouble()
    );
    Eigen::Matrix<double, 3, 1> velocity = Eigen::Matrix<double, 3, 1>(
        initialConditions["v0_sc_eci_m"][0].asDouble(),
        initialConditions["v0_sc_eci_m"][1].asDouble(),
        initialConditions["v0_sc_eci_m"][2].asDouble()
    );
    if (!initialConditions["t_J2000_s"].isNumeric()) {
        std::cerr << "simulation.initial_conditions.t_J2000_s must be seconds since J2000\n";
        return 1;
    }
    const double epoch_J2000_s = initialConditions["t_J2000_s"].asDouble();
    if (!std::isfinite(epoch_J2000_s)) {
        std::cerr << "simulation.initial_conditions.t_J2000_s must be finite\n";
        return 1;
    }
    double mass_kg = config["spacecraft"]["mass_kg"].asDouble();


    const Json::Value& gravityConfig = config["environment"]["gravity"];
    if (!gravityConfig.isObject() || !gravityConfig["model"].isString()) {
        std::cerr << "environment.gravity.model must be 'point_mass' or 'egm2008'\n";
        return 1;
    }

    const std::string gravityModelName = gravityConfig["model"].asString();
    GravityModel gravityModel;
    std::uint32_t gravityDegree = 0;
    std::uint32_t gravityOrder = 0;
    std::unique_ptr<EGM2008Coefficients> egm2008;
    if (gravityModelName == "point_mass") {
        gravityModel = GravityModel::E_PointMass;
    } else if (gravityModelName == "egm2008") {
        gravityModel = GravityModel::E_EGM2008;
        if (!gravityConfig["degree"].isUInt() || !gravityConfig["order"].isUInt()) {
            std::cerr << "environment.gravity.degree and order must be non-negative integers\n";
            return 1;
        }
        gravityDegree = gravityConfig["degree"].asUInt();
        gravityOrder = gravityConfig["order"].asUInt();
        try {
            egm2008 = std::make_unique<EGM2008Coefficients>("data/earth_gravity/emg2008_coefficients.bin");
        } catch (const std::exception& error) {
            std::cerr << "Unable to load EGM2008 coefficients: " << error.what() << '\n';
            return 1;
        }
        if (gravityDegree > egm2008->maxDegree || gravityOrder > gravityDegree) {
            std::cerr << "EGM2008 degree must not exceed " << egm2008->maxDegree
                      << " and order must not exceed degree\n";
            return 1;
        }
    } else {
        std::cerr << "Unsupported environment.gravity.model: " << gravityModelName << '\n';
        return 1;
    }
    Gravity gravity(egm2008.get());

    std::ofstream outputFile("output.csv");
    if (!outputFile) {
        std::cerr << "Unable to open output.csv for writing\n";
        return 1;
    }
    outputFile << "t_sec,r_x_m,r_y_m,r_z_m,v_x_m_s,v_y_m_s,v_z_m_s\n";
    
    // while loop for simulation
    double t_simulation {0.0};
    double sim_length_sec = config["simulation"]["sim_length_sec"].asDouble();
    double dt_sec = config["simulation"]["dt_sec"].asDouble();
    if (sim_length_sec <= 0.0 || dt_sec <= 0.0 || mass_kg <= 0.0) {
        std::cerr << "Simulation length, timestep, and spacecraft mass must be valid positive values\n";
        return 1;
    }

    while (B_simulationRunning) {
        const double t_J2000_s = epoch_J2000_s + t_simulation;
        const aocs::environment::JulianDate date =
            aocs::environment::JulianDate::fromJ2000Seconds(t_J2000_s);
        const auto earthRotation =
            aocs::environment::EarthRotation(date, date).reduce();
        const Eigen::Vector3d position_ecef = earthRotation.gcrsToItrs * position;
        const Eigen::Vector3d gravityForce_ecef =
            gravity.computeGravity(
                position_ecef, mass_kg, gravityModel, gravityDegree, gravityOrder);
        const Eigen::Vector3d gravityForce_eci =
            earthRotation.gcrsToItrs.transpose() * gravityForce_ecef;

        // update position and velocity
        velocity += gravityForce_eci / mass_kg * dt_sec;
        position += velocity * dt_sec;

        // update simulation time
        t_simulation += dt_sec;
        outputFile << t_simulation << ',' << position.x() << ',' << position.y() << ',' << position.z()
                   << ',' << velocity.x() << ',' << velocity.y() << ',' << velocity.z() << '\n';



        // check if simulation has reached the end time
        if (t_simulation >= sim_length_sec) {
            std::cout << "Simulation ended at time: " << t_simulation << " sec\n";
            std::cout << "Final position: " << position.transpose() << "\n";
            std::cout << "Final velocity: " << velocity.transpose() << "\n";
            B_simulationRunning = false;
        }

    }

    // clean up simulation resources
    outputFile.close();

    return 0;
}