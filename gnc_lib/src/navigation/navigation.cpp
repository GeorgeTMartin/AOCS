#include "navigation/navigation.hpp"

namespace aocs::navigation {

Navigation::Navigation(
    OrbitDeterminationMethod method,
    const Eigen::Vector3d&
        accelerationProcessNoiseSpectralDensityRTN_m2_s3)
    : orbitDetermination(
          method,
          accelerationProcessNoiseSpectralDensityRTN_m2_s3) {}

std::optional<InertialStateEstimate> Navigation::runStep(
    const NavigationInput& input) {
    if (input.gpsMeasurement) {
        if (hasInitialEstimate) {
            orbitDetermination.UpdateInertialState(
                *input.gpsMeasurement,
                input.gpsMeasurement->t_estimate_J2000_s);
        } else {
            orbitDetermination.InitializeInertialState(*input.gpsMeasurement);
            hasInitialEstimate = true;
        }
    }

    if (!hasInitialEstimate) {
        return std::nullopt;
    }

    orbitDetermination.PropagateTo(input.time_J2000_s);
    return orbitDetermination.GetInertialStateEstimate();
}

} // namespace aocs::navigation
