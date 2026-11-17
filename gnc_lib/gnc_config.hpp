#pragma once

#include "navigation/orbit_determination.hpp"
#include <Eigen/Dense>

namespace aocs {

class GncConfig {
public:
    static constexpr navigation::OrbitDeterminationMethod orbitDeterminationMethod =
        navigation::OrbitDeterminationMethod::ExtendedKalmanFilter;
    static inline const Eigen::Vector3d accelerationProcessNoiseSpectralDensityRTN_m2_s3 = Eigen::Vector3d(0.0001, 0.01, 0.0001);
};

} // namespace aocs