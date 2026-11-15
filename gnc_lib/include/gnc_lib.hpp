#pragma once

#include "ipc.hpp"
#include "navigation/orbit_determination.hpp"
#include "gnc_config.hpp"

namespace aocs {

struct GncState {
    GncState()
        : orbitDetermination(
              GncConfig::orbitDeterminationMethod,
              GncConfig::accelerationProcessNoiseSpectralDensity_m2_s3) {}

    OrbitDetermination orbitDetermination;
    bool hasInitialEstimate = false;
};

GncState initializeGnc();
int runGncStep(GncState& state, sharedMemoryMap& sharedMemory);

} // namespace aocs
