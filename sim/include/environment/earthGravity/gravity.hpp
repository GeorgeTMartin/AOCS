#pragma once

#include <Eigen/Dense>
#include <cassert>

class Gravity {
    public:
        Gravity();
        ~Gravity();
        Eigen::Matrix<double, 3, 1> computeGravity(const Eigen::Matrix<double, 3, 1>& position, double& mass_kg) const;
    private:
        double GM; // gravitational parameter
        double R;  // radius of the Earth
};

Gravity::Gravity() {
    GM = 3.986004418e14; // m^3/s^2
    R = 6.371e6;         // m
}

Gravity::~Gravity() {
}

Eigen::Matrix<double, 3, 1> Gravity::computeGravity(const Eigen::Matrix<double, 3, 1>& position, double& mass_kg) const {
    double r = position.norm();
    assert (r >= R); // prevent division by zero or negative values inside the Earth
    Eigen::Matrix<double, 3, 1> gravity;
    gravity = -GM / (r * r * r) * position * mass_kg;
    return gravity;
}