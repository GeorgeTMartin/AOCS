#include "ipc.hpp"
#include "environment/earthGravity/gravity.hpp"
#include <Eigen/Dense>
#include <json/value.h>
#include <json/reader.h>
#include <fstream>
#include <iostream>

bool B_simulationRunning = true;

int main() {
    // initialize simulation
    // load initial conditions from config.json
    std::ifstream configFile("config.json");
    if (!configFile) {
        std::cerr << "Unable to open config.json\n";
        return 1;
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
    double mass_kg = config["spacecraft"]["mass_kg"].asDouble();
    Gravity gravity;
    // while loop for simulation
    double t_simulation {0.0};
    double sim_length_sec = config["simulation"]["sim_length_sec"].asDouble();
    double dt_sec = config["simulation"]["dt_sec"].asDouble();
    if (sim_length_sec <= 0.0 || dt_sec <= 0.0 || mass_kg <= 0.0) {
        std::cerr << "Simulation length, timestep, and spacecraft mass must be valid positive values\n";
        return 1;
    }

    while (B_simulationRunning) {
        // update simulation state
        Eigen::Matrix<double, 3, 1> gravityForce = gravity.computeGravity(position, mass_kg);

        // update position and velocity
        velocity += gravityForce / mass_kg * dt_sec;
        position += velocity * dt_sec;

        // update simulation time
        t_simulation += dt_sec;

        // check if simulation has reached the end time
        if (t_simulation >= sim_length_sec) {
            std::cout << "Simulation ended at time: " << t_simulation << " sec\n";
            std::cout << "Final position: " << position.transpose() << "\n";
            std::cout << "Final velocity: " << velocity.transpose() << "\n";
            B_simulationRunning = false;
        }

    }

    // clean up simulation resources

    return 0;
}