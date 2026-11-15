#include "navigation/orbit_determination.hpp"

OrbitDetermination::OrbitDetermination(OrbitDeterminationMethod method)
    : method(method),
      inertialState(Eigen::VectorXd::Zero(6)),
      covariance(Eigen::MatrixXd::Identity(6, 6)) {
    if (method == OrbitDeterminationMethod::ExtendedKalmanFilter) {
        // Initialize EKF specific variables
    } else if (method == OrbitDeterminationMethod::UnscentedKalmanFilter) {
        // Initialize UKF specific variables
    }
}

void OrbitDetermination::UpdateInertialState(const Eigen::VectorXd& measurement, double& t_J2000_s) {
    // Implementation for updating the inertial state based on the measurement and time
}

const Eigen::VectorXd& OrbitDetermination::GetInertialState() const {
    return inertialState;
}

const Eigen::MatrixXd& OrbitDetermination::GetCovariance() const {
    return covariance;
}