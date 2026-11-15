#include "gnc_lib.hpp"
#include "gnc_config.hpp"
#include "environment/earth_rotation/iau2000b.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace aocs {
namespace {

constexpr std::size_t kStateDimension = ipc::kStateDimension;
using StateCovarianceRowMajor =
    Eigen::Matrix<double, kStateDimension, kStateDimension, Eigen::RowMajor>;

void publishNavigationEstimate(
    const InertialStateEstimate& estimate,
    navigationEstimateMessage& message) {
    message.valid = 0;
    if (estimate.covariance.rows() != static_cast<Eigen::Index>(kStateDimension) ||
        estimate.covariance.cols() != static_cast<Eigen::Index>(kStateDimension) ||
        !std::isfinite(estimate.t_estimate_J2000_s) ||
        !estimate.vec3_r_sc_eci.allFinite() ||
        !estimate.vec3_v_sc_eci.allFinite() ||
        !estimate.covariance.allFinite()) {
        return;
    }

    message.time_J2000_s = estimate.t_estimate_J2000_s;
    for (Eigen::Index axis = 0; axis < 3; ++axis) {
        message.position_ECI_m[static_cast<std::size_t>(axis)] =
            estimate.vec3_r_sc_eci(axis);
        message.velocity_ECI_m_s[static_cast<std::size_t>(axis)] =
            estimate.vec3_v_sc_eci(axis);
    }

    const StateCovarianceRowMajor covariance = estimate.covariance;
    std::copy(
        covariance.data(),
        covariance.data() + covariance.size(),
        message.covariance_row_major.begin());
    message.valid = 1;
}

} // namespace

GncState initializeGnc() {
    return {
        OrbitDetermination(GncConfig::orbitDeterminationMethod)
    };
}

int runGncStep(GncState& state, sharedMemoryMap& sharedMemory) {
    sharedMemory.navigationEstimate.valid = 0;

    const sensorMessage& sensors = sharedMemory.sensorData;
    const double t_J2000_s =
        sensors.aocs_time_J2000_s;
    double gpsTime_J2000_s =
        sensors.GPS_time_s - ipc::kGpsSecondsAtJ2000;
    const environment::JulianDate date = environment::JulianDate::fromJ2000Seconds(
        sensors.GPS_time_s - ipc::kGpsSecondsAtJ2000);
    const auto earthRotation = environment::EarthRotation(date, date).reduce();
    const Eigen::Matrix3d& gcrsToItrs = earthRotation.gcrsToItrs;
    const Eigen::Matrix3d& itrsToGcrs = gcrsToItrs.transpose();
    const Eigen::Vector3d earthAngularVelocity_ecef_rad_s(
        0.0, 0.0, 7.2921150e-5);

    InertialStateEstimate GPSMeasurement;
    GPSMeasurement.t_estimate_J2000_s = gpsTime_J2000_s;
    GPSMeasurement.vec3_r_sc_eci =
        itrsToGcrs * Eigen::Map<const Eigen::Vector3d>(
                          sensors.GPS_position_ECEF_m.data());
    const Eigen::Vector3d velocity_ecef_m_s =
        Eigen::Map<const Eigen::Vector3d>(sensors.GPS_velocity_ECEF_m_s.data());
    const Eigen::Vector3d position_ecef_m =
        Eigen::Map<const Eigen::Vector3d>(sensors.GPS_position_ECEF_m.data());
    GPSMeasurement.vec3_v_sc_eci =
        itrsToGcrs *
        (velocity_ecef_m_s + earthAngularVelocity_ecef_rad_s.cross(position_ecef_m));
    const Eigen::Map<const StateCovarianceRowMajor> measurementCovariance(
        sensors.GPS_covariance_row_major.data());
    GPSMeasurement.covariance = measurementCovariance;

    if (sensors.has_new_GPS_measurement != 0) {
        if (state.hasInitialEstimate) {
            state.orbitDetermination.UpdateInertialState(
                GPSMeasurement, gpsTime_J2000_s);
        } else {
            state.orbitDetermination.InitializeInertialState(GPSMeasurement);
            state.hasInitialEstimate = true;
        }
    } else {
        if (!state.hasInitialEstimate) {
            return 1;
        }
    }
    if (!state.hasInitialEstimate) {
        return 1;
    }
    state.orbitDetermination.PropagateTo(t_J2000_s);
    publishNavigationEstimate(
        state.orbitDetermination.GetInertialStateEstimate(),
        sharedMemory.navigationEstimate);

    return sharedMemory.navigationEstimate.valid == 1 ? 0 : 1;
}

} // namespace aocs