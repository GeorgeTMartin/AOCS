#pragma once

#include "ipc.hpp"
#include "gnc_config.hpp"
#include "navigation/navigation.hpp"

namespace aocs {

enum class GNCState {
    AwaitingNavigationInitialization,
    Operational
};

class GNCStateMachine {
public:
    void evaluate(bool hasNewGpsMeasurement) {
        if (state == GNCState::AwaitingNavigationInitialization &&
            hasNewGpsMeasurement) {
            state = GNCState::Operational;
        }
    }

    bool navigationReady() const {
        return state == GNCState::Operational;
    }

private:
    GNCState state = GNCState::AwaitingNavigationInitialization;
};

struct GNCExecutive {
    GNCExecutive()
        : navigation(
              GncConfig::orbitDeterminationMethod,
              GncConfig::accelerationProcessNoiseSpectralDensityRTN_m2_s3) {}

    GNCStateMachine stateMachine;
    navigation::Navigation navigation;
};

GNCExecutive initializeGnc();
int runGncStep(GNCExecutive& executive, sharedMemoryMap& sharedMemory);

} // namespace aocs
