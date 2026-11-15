#pragma once

#include <Eigen/Dense>

enum class OrbitDeterminationMethod {
    ExtendedKalmanFilter,
    UnscentedKalmanFilter
};

struct InertialStateEstimate {
  double t_estimate_J2000_s;
  Eigen::Vector3d vec3_r_sc_eci;
  Eigen::Vector3d vec3_v_sc_eci;
  Eigen::MatrixXd covariance;
};

class OrbitDetermination {
public:
    explicit OrbitDetermination(
        OrbitDeterminationMethod method,
        double accelerationProcessNoiseSpectralDensity_m2_s3);
    void InitializeInertialState(const InertialStateEstimate& initialEstimate);
    void PropagateTo(double t_J2000_s);
    void UpdateInertialState(InertialStateEstimate& measurement, double& t_J2000_s);
    const InertialStateEstimate& GetInertialStateEstimate() const;
    Eigen::Vector3d& calculateAcceleration(const Eigen::Vector3d& vec3_r_sc_eci, const Eigen::Vector3d& vec3_v_sc_eci, double t_start_J2000_s, double t_end_J2000_s);

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
    InertialStateEstimate measurement;
};