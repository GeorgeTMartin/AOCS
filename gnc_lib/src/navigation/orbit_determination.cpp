#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

#include "navigation/orbit_determination.hpp"
#include "flight_constants.hpp"

namespace aocs::navigation {

OrbitDetermination::OrbitDetermination(
    OrbitDeterminationMethod method,
    const Eigen::Vector3d& accelerationProcessNoiseSpectralDensityRTN_m2_s3)
    : method(method),
      accelerationProcessNoiseSpectralDensityRTN_m2_s3(accelerationProcessNoiseSpectralDensityRTN_m2_s3),
      inertialState{
          0.0,
          Eigen::Matrix<double,3,1>::Zero(),
          Eigen::Matrix<double,3,1>::Zero(),
          StateCovariance::Identity()} {
    if (!accelerationProcessNoiseSpectralDensityRTN_m2_s3.allFinite() ||
        accelerationProcessNoiseSpectralDensityRTN_m2_s3.minCoeff() < 0.0) {
        throw std::invalid_argument(
            "Acceleration process-noise spectral density must be finite and non-negative");
    }
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

void OrbitDetermination::PropagateInertialState(
    InertialStateEstimate& state,
    double t_start_J2000_s,
    double t_end_J2000_s) {
    const Eigen::Matrix<double,3,1> positionAtStart_eci_m = state.vec3_r_sc_eci;
    const Eigen::Matrix<double,3,1> velocityAtStart_eci_m_s = state.vec3_v_sc_eci;
    const StateCovariance covarianceAtStart = state.covariance;
    const double radius_m = positionAtStart_eci_m.norm();
    const double elapsedTime_s = t_end_J2000_s - t_start_J2000_s;
    if (!positionAtStart_eci_m.allFinite() ||
        !velocityAtStart_eci_m_s.allFinite() ||
        !std::isfinite(radius_m) || radius_m <= 0.0) {
        throw std::invalid_argument(
            "RTN process noise requires a finite, nonzero position");
    }
    const Eigen::Vector3d angularMomentum_eci =
        positionAtStart_eci_m.cross(velocityAtStart_eci_m_s);
    const double angularMomentumMagnitude =
        angularMomentum_eci.norm();
    if (!std::isfinite(angularMomentumMagnitude) ||
        angularMomentumMagnitude <= 0.0) {
        throw std::invalid_argument(
            "RTN process noise requires nonzero orbital angular momentum");
    }

    const Eigen::Matrix<double,3,1> vec3_a_sc_eci = -flight_constants::kEarthMu * positionAtStart_eci_m / (radius_m * radius_m * radius_m);

    state.vec3_r_sc_eci = positionAtStart_eci_m + velocityAtStart_eci_m_s * elapsedTime_s;
    state.vec3_v_sc_eci = velocityAtStart_eci_m_s + vec3_a_sc_eci * elapsedTime_s;
    state.t_estimate_J2000_s = t_end_J2000_s;

    const double radiusSquared_m2 = radius_m * radius_m;
    const double radiusCubed_m3 = radiusSquared_m2 * radius_m;
    const Eigen::Matrix3d accelerationGradient = flight_constants::kEarthMu * (3.0 * positionAtStart_eci_m * positionAtStart_eci_m.transpose() /
                                                         (radiusCubed_m3 * radiusSquared_m2) - Eigen::Matrix3d::Identity() / radiusCubed_m3);
    Eigen::Matrix<double, 6, 6> dynamicsJacobian = Eigen::Matrix<double, 6, 6>::Zero();
    dynamicsJacobian.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    dynamicsJacobian.block<3, 3>(3, 0) = accelerationGradient;

    Eigen::Matrix<double, 6, 6> stateTransitionMatrix =
        Eigen::Matrix<double, 6, 6>::Identity();
    stateTransitionMatrix += dynamicsJacobian * elapsedTime_s;
    const Eigen::Vector3d radial_eci = positionAtStart_eci_m / radius_m;
    const Eigen::Vector3d normal_eci =
        angularMomentum_eci / angularMomentumMagnitude;
    Eigen::Matrix3d rtnToEci;
    rtnToEci.col(0) = radial_eci;
    rtnToEci.col(1) = normal_eci.cross(radial_eci);
    rtnToEci.col(2) = normal_eci;
    const Eigen::Matrix3d processNoiseECI =
        rtnToEci *
        accelerationProcessNoiseSpectralDensityRTN_m2_s3.asDiagonal() *
        rtnToEci.transpose();

    Eigen::Matrix<double, 6, 6> processNoise =
        Eigen::Matrix<double, 6, 6>::Zero();
    processNoise.block<3, 3>(0, 0) = processNoiseECI * (elapsedTime_s * elapsedTime_s * elapsedTime_s / 3.0);
    processNoise.block<3, 3>(0, 3) = processNoiseECI * (elapsedTime_s * elapsedTime_s / 2.0);
    processNoise.block<3, 3>(3, 0) = processNoiseECI * (elapsedTime_s * elapsedTime_s / 2.0);
    processNoise.block<3, 3>(3, 3) = processNoiseECI * elapsedTime_s;

    state.covariance = stateTransitionMatrix * covarianceAtStart * stateTransitionMatrix.transpose() + processNoise;

    return;
}

void OrbitDetermination::MeasurementUpdate(
    InertialStateEstimate& estimate,
    const InertialStateEstimate& measurement) {
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

void OrbitDetermination::UpdateInertialState(
    const InertialStateEstimate& measurement,
    double t_J2000_s) {
    if (!std::isfinite(measurement.t_estimate_J2000_s) ||
        !std::isfinite(t_J2000_s) ||
        t_J2000_s < measurement.t_estimate_J2000_s) {
        throw std::invalid_argument(
            "Orbit determination update time must be finite and not precede the measurement");
    }

    PropagateInertialState(
        inertialState, inertialState.t_estimate_J2000_s, t_J2000_s);
    InertialStateEstimate propagatedMeasurement = measurement;
    PropagateInertialState(
        propagatedMeasurement,
        measurement.t_estimate_J2000_s,
        t_J2000_s);
    MeasurementUpdate(inertialState, propagatedMeasurement);
}

const InertialStateEstimate& OrbitDetermination::GetInertialStateEstimate() const {
    return inertialState;
}

} // namespace aocs::navigation