#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

#include "navigation/orbit_determination.hpp"
#include "flight_constants.hpp"


OrbitDetermination::OrbitDetermination(OrbitDeterminationMethod method)
    : method(method),
      inertialState{0.0, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::MatrixXd::Identity(6, 6)},
      measurement{0.0, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::MatrixXd::Identity(6, 6)} {
    if (method == OrbitDeterminationMethod::ExtendedKalmanFilter) {
        // Initialize EKF specific variables, if any
    } else if (method == OrbitDeterminationMethod::UnscentedKalmanFilter) {
        // Initialize UKF specific variables, if any
    }
}

void OrbitDetermination::InitializeInertialState(
    const InertialStateEstimate& initialEstimate) {
    inertialState = initialEstimate;
}

void OrbitDetermination::PropagateTo(double t_J2000_s) {
    if (!std::isfinite(t_J2000_s) ||
        t_J2000_s < inertialState.t_estimate_J2000_s) {
        throw std::invalid_argument(
            "Orbit estimate propagation time must be finite and nondecreasing");
    }
    PropagateInertialState(
        inertialState, inertialState.t_estimate_J2000_s, t_J2000_s);
}

void OrbitDetermination::PropagateInertialState(const InertialStateEstimate& state, double t_start_J2000_s, double t_end_J2000_s) {
    // Propagate state with nonlinear dynamics

    const double radius_m = state.vec3_r_sc_eci.norm();
    const Eigen::Vector3d vec3_a_sc_eci =
        -flight_constants::kEarthMu * state.vec3_r_sc_eci /
        (radius_m * radius_m * radius_m);

    inertialState.vec3_r_sc_eci += state.vec3_v_sc_eci * (t_end_J2000_s - t_start_J2000_s);
    inertialState.vec3_v_sc_eci += vec3_a_sc_eci * (t_end_J2000_s - t_start_J2000_s);
    inertialState.t_estimate_J2000_s = t_end_J2000_s;

    // Linearize and Propagate Covariance
    Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(6, 6);
    Jacobian.block(0,3,3,3) = Eigen::Matrix3d::Identity();
    Jacobian(3,0) = flight_constants::kEarthMu * (3.0 * state.vec3_r_sc_eci.x() * state.vec3_r_sc_eci.x() / state.vec3_r_sc_eci.squaredNorm() - 1.0) / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(3,1) = flight_constants::kEarthMu * 3.0 * state.vec3_r_sc_eci.x() * state.vec3_r_sc_eci.y() / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(3,2) = flight_constants::kEarthMu * 3.0 * state.vec3_r_sc_eci.x() * state.vec3_r_sc_eci.z() / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(4,0) = flight_constants::kEarthMu * 3.0 * state.vec3_r_sc_eci.x() * state.vec3_r_sc_eci.y() / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(4,1) = flight_constants::kEarthMu * (3.0 * state.vec3_r_sc_eci.y() * state.vec3_r_sc_eci.y() / state.vec3_r_sc_eci.squaredNorm() - 1.0) / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(4,2) = flight_constants::kEarthMu * 3.0 * state.vec3_r_sc_eci.y() * state.vec3_r_sc_eci.z() / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(5,0) = flight_constants::kEarthMu * 3.0 * state.vec3_r_sc_eci.x() * state.vec3_r_sc_eci.z() / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(5,1) = flight_constants::kEarthMu * 3.0 * state.vec3_r_sc_eci.y() * state.vec3_r_sc_eci.z() / state.vec3_r_sc_eci.squaredNorm();
    Jacobian(5,2) = flight_constants::kEarthMu * (3.0 * state.vec3_r_sc_eci.z() * state.vec3_r_sc_eci.z() / state.vec3_r_sc_eci.squaredNorm() - 1.0) / state.vec3_r_sc_eci.squaredNorm();

    Eigen::MatrixXd StateTransitionMatrix = Eigen::MatrixXd::Identity(6, 6);
    StateTransitionMatrix += Jacobian * (t_end_J2000_s - t_start_J2000_s);
    inertialState.covariance = StateTransitionMatrix * state.covariance * StateTransitionMatrix.transpose();

    return;
}

void OrbitDetermination::MeasurementUpdate(InertialStateEstimate& estimate, const InertialStateEstimate& measurement) {
    Eigen::Matrix<double, 6, 1> estimatedState;
    estimatedState << estimate.vec3_r_sc_eci, estimate.vec3_v_sc_eci;
    Eigen::Matrix<double, 6, 1> measuredState;
    measuredState << measurement.vec3_r_sc_eci, measurement.vec3_v_sc_eci;

    const Eigen::Matrix<double, 6, 6> innovationCovariance =
        estimate.covariance + measurement.covariance;
    const Eigen::Matrix<double, 6, 6> kalmanGain =
        estimate.covariance * innovationCovariance.inverse();
    const Eigen::Matrix<double, 6, 1> updatedState =
        estimatedState + kalmanGain * (measuredState - estimatedState);

    estimate.vec3_r_sc_eci = updatedState.head<3>();
    estimate.vec3_v_sc_eci = updatedState.tail<3>();
    estimate.covariance =
        (Eigen::Matrix<double, 6, 6>::Identity() - kalmanGain) * estimate.covariance;
}

void OrbitDetermination::UpdateInertialState(InertialStateEstimate& GPSMeasurement, double& t_J2000_s) {
    // Implementation for updating the inertial state based on the measurement and time
    // Propagate the inertial state and covariance to the new time
    PropagateInertialState(inertialState, inertialState.t_estimate_J2000_s, t_J2000_s);
    // Propagate the measurement and measurement covariance to the new time
    PropagateInertialState(GPSMeasurement, t_J2000_s, inertialState.t_estimate_J2000_s);
    // measurement update
    MeasurementUpdate(inertialState, GPSMeasurement);
}

const InertialStateEstimate& OrbitDetermination::GetInertialStateEstimate() const {
    return inertialState;
}