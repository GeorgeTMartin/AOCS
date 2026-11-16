#pragma once

#include <Eigen/Dense>

namespace aocs::navigation {

enum class OrbitDeterminationMethod {
    ExtendedKalmanFilter,
    UnscentedKalmanFilter
};

struct InertialStateEstimate {
    double t_estimate_J2000_s = 0.0;
    Eigen::Vector3d vec3_r_sc_eci = Eigen::Vector3d::Zero();
    Eigen::Vector3d vec3_v_sc_eci = Eigen::Vector3d::Zero();
    Eigen::MatrixXd covariance;
};

class OrbitDetermination {
public:
    explicit OrbitDetermination(
        OrbitDeterminationMethod method,
        double accelerationProcessNoiseSpectralDensity_m2_s3);

    void InitializeInertialState(const InertialStateEstimate& initialEstimate);
    void PropagateTo(double t_J2000_s);
    void UpdateInertialState(
        const InertialStateEstimate& measurement,
        double t_J2000_s);
    const InertialStateEstimate& GetInertialStateEstimate() const;

private:
    void PropagateInertialState(
        InertialStateEstimate& state,
        double t_start_J2000_s,
        double t_end_J2000_s);
    void MeasurementUpdate(
        InertialStateEstimate& estimate,
        const InertialStateEstimate& measurement);

    OrbitDeterminationMethod method;
    double accelerationProcessNoiseSpectralDensity_m2_s3;
    InertialStateEstimate inertialState;
};

} // namespace aocs::navigation
