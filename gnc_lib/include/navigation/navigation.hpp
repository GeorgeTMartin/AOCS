#pragma once

#include "navigation/orbit_determination.hpp"
#include <optional>

namespace aocs::navigation {

struct NavigationInput {
    double time_J2000_s = 0.0;
    std::optional<InertialStateEstimate> gpsMeasurement;
};

class Navigation {
public:
    Navigation(
        OrbitDeterminationMethod method,
        double accelerationProcessNoiseSpectralDensity_m2_s3);

    std::optional<InertialStateEstimate> runStep(
        const NavigationInput& input);

private:
    OrbitDetermination orbitDetermination;
    bool hasInitialEstimate = false;
};

} // namespace aocs::navigation
