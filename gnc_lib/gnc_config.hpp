#pragma once

#include "navigation/orbit_determination.hpp"

namespace aocs {

class GncConfig {
public:
    static constexpr navigation::OrbitDeterminationMethod orbitDeterminationMethod =
        navigation::OrbitDeterminationMethod::ExtendedKalmanFilter;
    static constexpr double accelerationProcessNoiseSpectralDensity_m2_s3 = 1.0e-4;
};

} // namespace aocs