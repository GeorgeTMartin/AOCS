#pragma once

#include "include/navigation/orbit_determination.hpp"

class GncConfig {
public:
    static inline OrbitDeterminationMethod orbitDeterminationMethod = OrbitDeterminationMethod::ExtendedKalmanFilter;
};