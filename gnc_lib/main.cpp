#include "gnc_lib.hpp"
#include "conversions/frames/FK5.hpp"
#include "flight_constants.hpp"
#include "navigation/navigation.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace aocs {
namespace {

navigation::NavigationInput mapNavigationInput(const sensorMessage& sensors) {
    navigation::NavigationInput input;
    input.time_J2000_s = sensors.aocs_time_J2000_s;
    if (sensors.has_new_GPS_measurement == 0) {
        return input;
    }

    double gpsTime_J2000_s =
        sensors.GPS_time_s - flight_constants::kGPSSecondsAtJ2000;
    const double gpsTimeResolution_s =
        4.0 * std::numeric_limits<double>::epsilon() *
        std::max(std::abs(sensors.GPS_time_s),
                 flight_constants::kGPSSecondsAtJ2000);
    if (std::abs(gpsTime_J2000_s - input.time_J2000_s) <=
        gpsTimeResolution_s) {
        gpsTime_J2000_s = input.time_J2000_s;
    }

    const Eigen::Quaterniond earthRotation =
        getQuaternionECEFtoECI(gpsTime_J2000_s);
    const Eigen::Matrix3d ecefToEci = earthRotation.toRotationMatrix();
    const Eigen::Vector3d position_ecef_m = Eigen::Map<const Eigen::Vector3d>(sensors.GPS_position_ECEF_m.data());
    const Eigen::Vector3d velocity_ecef_m_s = Eigen::Map<const Eigen::Vector3d>(sensors.GPS_velocity_ECEF_m_s.data());
    const Eigen::Vector3d earthAngularVelocity_ecef_rad_s(0.0, 0.0, flight_constants::kEarthAngularRate);

    navigation::InertialStateEstimate measurement;
    measurement.t_estimate_J2000_s = gpsTime_J2000_s;
    measurement.vec3_r_sc_eci = ecefToEci * position_ecef_m;
    measurement.vec3_v_sc_eci = ecefToEci * (velocity_ecef_m_s + earthAngularVelocity_ecef_rad_s.cross(position_ecef_m));

    measurement.covariance =
        Eigen::Map<const navigation::StateCovariance>(
            sensors.GPS_covariance_row_major.data());
    input.gpsMeasurement = std::move(measurement);
    return input;
}

void mapNavigationOutputToIpc(
    const navigation::InertialStateEstimate& estimate,
    navigationEstimateMessage& message) {
    message.valid = 0;
    if (!std::isfinite(estimate.t_estimate_J2000_s) ||
        !estimate.vec3_r_sc_eci.allFinite() ||
        !estimate.vec3_v_sc_eci.allFinite() ||
        !estimate.covariance.allFinite()) {
        return;
    }

    message.time_J2000_s = estimate.t_estimate_J2000_s;
    for (Eigen::Index axis = 0;
         axis < static_cast<Eigen::Index>(ipc::kStateDimension / 2);
         ++axis) {
        const std::size_t index = static_cast<std::size_t>(axis);
        message.position_ECI_m[index] = estimate.vec3_r_sc_eci(axis);
        message.velocity_ECI_m_s[index] = estimate.vec3_v_sc_eci(axis);
    }

    const navigation::StateCovariance covariance = estimate.covariance;
    std::copy(
        covariance.data(),
        covariance.data() + covariance.size(),
        message.covariance_row_major.begin());
    message.valid = 1;
}

} // namespace

GNCExecutive initializeGnc() {
    return GNCExecutive();
}

int runGncStep(GNCExecutive& executive, sharedMemoryMap& sharedMemory) {
    navigationEstimateMessage& navigationMessage =
        sharedMemory.navigationEstimate;
    navigationMessage.valid = 0;

    const navigation::NavigationInput inputs =
        mapNavigationInput(sharedMemory.sensorData);
    executive.stateMachine.evaluate(inputs.gpsMeasurement.has_value());
    if (!executive.stateMachine.navigationReady()) {
        return 1;
    }

    const std::optional<navigation::InertialStateEstimate> estimate =
        executive.navigation.runStep(inputs);
    if (!estimate) {
        return 1;
    }

    mapNavigationOutputToIpc(*estimate, navigationMessage);
    return navigationMessage.valid == 1 ? 0 : 1;
}

} // namespace aocs
