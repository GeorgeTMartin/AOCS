#pragma once

#include "ipc.hpp"
#include "navigation/orbit_determination.hpp"
#include "gnc_config.hpp"

namespace aocs {

struct GncState {
    explicit GncState(
        double accelerationProcessNoiseSpectralDensity_m2_s3 =
            GncConfig::accelerationProcessNoiseSpectralDensity_m2_s3)
        : orbitDetermination(
              GncConfig::orbitDeterminationMethod,
              accelerationProcessNoiseSpectralDensity_m2_s3) {}

    OrbitDetermination orbitDetermination;
    bool hasInitialEstimate = false;
};

GncState initializeGnc(
    double accelerationProcessNoiseSpectralDensity_m2_s3 =
        GncConfig::accelerationProcessNoiseSpectralDensity_m2_s3);
int runGncStep(GncState& state, sharedMemoryMap& sharedMemory);

} // namespace aocs
