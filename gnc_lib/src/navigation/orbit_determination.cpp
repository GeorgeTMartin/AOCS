#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

#include "navigation/orbit_determination.hpp"
#include "flight_constants.hpp"


OrbitDetermination::OrbitDetermination(
    OrbitDeterminationMethod method,
    double accelerationProcessNoiseSpectralDensity_m2_s3)
    : method(method),
      accelerationProcessNoiseSpectralDensity_m2_s3(
          accelerationProcessNoiseSpectralDensity_m2_s3),
      inertialState{0.0, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::MatrixXd::Identity(6, 6)},
      measurement{0.0, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::MatrixXd::Identity(6, 6)} {
    if (!std::isfinite(accelerationProcessNoiseSpectralDensity_m2_s3) ||
        accelerationProcessNoiseSpectralDensity_m2_s3 < 0.0) {
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
    const Eigen::Vector3d positionAtStart_eci_m = state.vec3_r_sc_eci;
    const Eigen::Vector3d velocityAtStart_eci_m_s = state.vec3_v_sc_eci;
    const Eigen::MatrixXd covarianceAtStart = state.covariance;
    const double radius_m = positionAtStart_eci_m.norm();
    const double elapsedTime_s = t_end_J2000_s - t_start_J2000_s;
    const Eigen::Vector3d vec3_a_sc_eci =
        -flight_constants::kEarthMu * positionAtStart_eci_m /
        (radius_m * radius_m * radius_m);

    state.vec3_r_sc_eci =
        positionAtStart_eci_m + velocityAtStart_eci_m_s * elapsedTime_s;
    state.vec3_v_sc_eci =
        velocityAtStart_eci_m_s + vec3_a_sc_eci * elapsedTime_s;
    state.t_estimate_J2000_s = t_end_J2000_s;

    const double radiusSquared_m2 = radius_m * radius_m;
    const double radiusCubed_m3 = radiusSquared_m2 * radius_m;
    const Eigen::Matrix3d accelerationGradient =
        flight_constants::kEarthMu *
        (3.0 * positionAtStart_eci_m * positionAtStart_eci_m.transpose() /
             (radiusCubed_m3 * radiusSquared_m2) -
         Eigen::Matrix3d::Identity() / radiusCubed_m3);
    Eigen::Matrix<double, 6, 6> dynamicsJacobian =
        Eigen::Matrix<double, 6, 6>::Zero();
    dynamicsJacobian.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
    dynamicsJacobian.block<3, 3>(3, 0) = accelerationGradient;

    Eigen::Matrix<double, 6, 6> stateTransitionMatrix =
        Eigen::Matrix<double, 6, 6>::Identity();
    stateTransitionMatrix += dynamicsJacobian * elapsedTime_s;
    Eigen::Matrix<double, 6, 6> processNoise =
        Eigen::Matrix<double, 6, 6>::Zero();
    const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
    const double spectralDensity =
        accelerationProcessNoiseSpectralDensity_m2_s3;
    processNoise.block<3, 3>(0, 0) =
        identity * (spectralDensity * elapsedTime_s * elapsedTime_s *
                    elapsedTime_s / 3.0);
    processNoise.block<3, 3>(0, 3) =
        identity * (spectralDensity * elapsedTime_s * elapsedTime_s / 2.0);
    processNoise.block<3, 3>(3, 0) =
        identity * (spectralDensity * elapsedTime_s * elapsedTime_s / 2.0);
    processNoise.block<3, 3>(3, 3) =
        identity * (spectralDensity * elapsedTime_s);
    state.covariance =
        stateTransitionMatrix * covarianceAtStart *
            stateTransitionMatrix.transpose() +
        processNoise;

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