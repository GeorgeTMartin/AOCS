#pragma once

#include "include/navigation/orbit_determination.hpp"

class GncConfig {
public:
    static inline OrbitDeterminationMethod orbitDeterminationMethod = OrbitDeterminationMethod::ExtendedKalmanFilter;
    static constexpr double accelerationProcessNoiseSpectralDensity_m2_s3 = 1.0e-4;
};