#pragma once

#include <Eigen/Dense>

enum class OrbitDeterminationMethod {
    ExtendedKalmanFilter,
    UnscentedKalmanFilter
};


class OrbitDetermination {
public:
    explicit OrbitDetermination(OrbitDeterminationMethod method);
    void UpdateInertialState(const Eigen::VectorXd& measurement, double& t_J2000_s);
    Eigen::VectorXd& GetInertialState();
    const Eigen::VectorXd& GetInertialState() const;
    Eigen::MatrixXd& GetCovariance();
    const Eigen::MatrixXd& GetCovariance() const;

private:
    OrbitDeterminationMethod method;
    Eigen::VectorXd inertialState;
    Eigen::MatrixXd covariance;
};