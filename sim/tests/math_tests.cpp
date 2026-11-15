#include "environment/earth_rotation/iau2000b.hpp"
#include "gnc_lib.hpp"
#include "ipc.hpp"
#include <math/quaternion.hpp>
#include "simulation/simulation.hpp"
#include "sim_constants.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void requireNear(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        std::ostringstream error;
        error << message << " (actual " << std::setprecision(17) << actual
              << ", expected " << expected << ')';
        throw std::runtime_error(error.str());
    }
}

void testQuaternion() {
    using Quaternion = sim::math::Quaternion<double>;
    using Vector3 = Eigen::Vector3d;

    static_assert(sizeof(Quaternion) == 4 * sizeof(double));

    const Quaternion q = Quaternion::fromPassiveAxisAngle(
        Vector3::UnitZ(), 0.5 * SimConstants::kPi);
    const Vector3 rotated = q.rotatePassive(Vector3::UnitX());
    require((rotated - Vector3(0.0, -1.0, 0.0)).norm() < 1e-14,
            "Passive +90 degree Z rotation must map X to -Y");
    require((q.toPassiveMatrix() * Vector3::UnitX() - rotated).norm() < 1e-14,
            "Quaternion vector and matrix rotation must agree");
    require((Quaternion::fromPassiveMatrix(q.toPassiveMatrix()).toPassiveMatrix() -
             q.toPassiveMatrix()).norm() < 1e-14,
            "Quaternion matrix conversion must round-trip");

    const Quaternion product = q * q;
    require((product.toPassiveMatrix() -
             q.toPassiveMatrix() * q.toPassiveMatrix()).norm() < 1e-14,
            "Hamilton quaternion composition must match matrix composition");
    const Quaternion qx =
        Quaternion::fromPassiveAxisAngle(Vector3::UnitX(), 0.3);
    const Quaternion qy =
        Quaternion::fromPassiveAxisAngle(Vector3::UnitY(), -0.7);
    require(((qx * qy).toPassiveMatrix() -
             qy.toPassiveMatrix() * qx.toPassiveMatrix()).norm() < 1e-14,
            "Hamilton composition must preserve passive rotation ordering");

    bool rejectedZero = false;
    try {
        (void)Quaternion(0.0, 0.0, 0.0, 0.0).normalized();
    } catch (const std::invalid_argument&) {
        rejectedZero = true;
    }
    require(rejectedZero, "Normalizing a zero quaternion must fail");
}

void testIau2000B() {
    using aocs::environment::EarthRotation;
    using aocs::environment::JulianDate;

    const JulianDate date{2400000.5, 53736.0};
    const auto reduction = EarthRotation(date, date).reduce();
    requireNear(reduction.nutation.longitude, -0.9632552291148362783e-5,
                1e-13, "IAU 2000B nutation longitude reference mismatch");
    requireNear(reduction.nutation.obliquity, 0.4063197106621159367e-4,
                1e-13, "IAU 2000B nutation obliquity reference mismatch");
    requireNear(reduction.greenwichMeanSiderealTime, 1.754174972210740592,
                1e-12, "IAU 2000 Greenwich mean sidereal time reference mismatch");
    requireNear(reduction.equationOfEquinoxes, -0.8835700060003032831e-5,
                1e-18, "IAU 2000B equation-of-equinoxes reference mismatch");

    const JulianDate matrixDate{2400000.5, 50123.9999};
    const auto matrixReduction = EarthRotation(matrixDate, matrixDate).reduce();
    const Eigen::Matrix3d& gcrsToTrueEquinox = matrixReduction.gcrsToTrueEquinox;
    requireNear(gcrsToTrueEquinox(0, 0), 0.9999995832776208280,
                1e-12, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(0, 1), 0.8372401264429654837e-3,
                1e-14, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(0, 2), 0.3639691681450271771e-3,
                1e-14, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(1, 0), -0.8372552234147137424e-3,
                1e-14, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(1, 1), 0.9999996486477686123,
                1e-12, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(1, 2), 0.4132832190946052890e-4,
                1e-14, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(2, 0), -0.3639344385341866407e-3,
                1e-14, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(2, 1), -0.4163303977421522785e-4,
                1e-14, "IAU 2000B bias-precession-nutation reference mismatch");
    requireNear(gcrsToTrueEquinox(2, 2), 0.9999999329092049734,
                1e-12, "IAU 2000B bias-precession-nutation reference mismatch");

    const Eigen::Matrix3d& matrix = reduction.gcrsToItrs;
    require((matrix.transpose() * matrix - Eigen::Matrix3d::Identity()).norm() < 1e-12,
            "ECI-to-ECEF matrix must be orthogonal");
    requireNear(matrix.determinant(), 1.0, 1e-12,
                "ECI-to-ECEF matrix must have unit determinant");
    require((reduction.gcrsToItrsQuaternion.toPassiveMatrix() - matrix).norm() < 1e-12,
            "ECI-to-ECEF quaternion and matrix must agree");

    const JulianDate tt = JulianDate::fromJ2000Seconds(86400.0);
    requireNear(tt.part1, 2451545.0, 0.0,
                "J2000-seconds date conversion must retain the split epoch");
    requireNear(tt.part2, 1.0, 0.0,
                "J2000-seconds date conversion must convert seconds to days");

    bool rejectedNonFinite = false;
    try {
        (void)EarthRotation({2451545.0, 0.0}, {2451545.0, 0.0},
                            std::numeric_limits<double>::quiet_NaN(), 0.0);
    } catch (const std::invalid_argument&) {
        rejectedNonFinite = true;
    }
    require(rejectedNonFinite, "Non-finite polar motion must be rejected");
}

void testRK4Propagation() {
    using aocs::simulation::Propagation;
    using aocs::simulation::SpacecraftState;

    SpacecraftState state{
        Eigen::Vector3d::UnitX(),
        Eigen::Vector3d::Zero()};
    double elapsedTime_s = 0.0;
    const double dt_sec = 0.1;
    Propagation propagation;
    propagation.propagate(
        state,
        elapsedTime_s,
        dt_sec,
        [](const SpacecraftState& intermediate, double) -> Eigen::Vector3d {
            return Eigen::Vector3d(-intermediate.position_eci_m);
        });

    requireNear(state.position_eci_m.x(), std::cos(dt_sec), 1e-7,
                "RK4 harmonic oscillator position mismatch");
    requireNear(state.velocity_eci_m_s.x(), -std::sin(dt_sec), 1e-7,
                "RK4 harmonic oscillator velocity mismatch");
    requireNear(elapsedTime_s, dt_sec, 0.0,
                "RK4 propagation must advance elapsed simulation time");
}

sensorMessage receivedSensorMessage;

int captureSensorMessage(sharedMemoryMap& memory) {
    receivedSensorMessage = memory.sensorData;
    return 0;
}

void testSharedMemoryGncTransport() {
    const std::string name =
        "/aocs_ipc_test_" + std::to_string(static_cast<long long>(getpid()));
    SharedMemoryRegion region(name);
    SharedMemoryGncTransport transport(region.data(), &captureSensorMessage);

    sensorMessage input;
    input.GPS_position_ECEF_m = {1.0, 2.0, 3.0};
    input.GPS_velocity_ECEF_m_s = {4.0, 5.0, 6.0};
    input.GPS_time_s = 7.0;
    input.aocs_time_J2000_s = 6.0;
    input.has_new_GPS_measurement = 1;
    input.GPS_covariance_row_major[0] = 9.0;
    require(transport.exchange(input) == 0,
            "Shared-memory GNC transport must return the GNC status");
    require(receivedSensorMessage.GPS_position_ECEF_m ==
                input.GPS_position_ECEF_m,
            "GNC entry point must receive the published GPS position");
    require(receivedSensorMessage.GPS_velocity_ECEF_m_s ==
                input.GPS_velocity_ECEF_m_s,
            "GNC entry point must receive the published GPS velocity");
    requireNear(receivedSensorMessage.GPS_time_s, input.GPS_time_s, 0.0,
                "GNC entry point must receive the published GPS time");
    requireNear(
        receivedSensorMessage.GPS_covariance_row_major[0],
        input.GPS_covariance_row_major[0],
        0.0,
        "GNC entry point must receive the GPS covariance");
}

void testSharedMemoryRegionRecoversStaleObject() {
    const std::string name =
        "/aocs_stale_shm_test_" + std::to_string(static_cast<long long>(getpid()));
    const int staleDescriptor = shm_open(
        name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
    require(staleDescriptor != -1,
            "Test must create a stale shared-memory object");
    require(ftruncate(staleDescriptor, sizeof(sharedMemoryMap)) == 0,
            "Test must size its stale shared-memory object");
    close(staleDescriptor);

    {
        SharedMemoryRegion recoveredRegion(name);
        require(recoveredRegion.data().navigationEstimate.valid == 0,
        "Stale shared-memory contents must be reinitialized");
        bool rejectedConcurrentOwner = false;
        try {
            SharedMemoryRegion concurrentRegion(name);
        } catch (const std::runtime_error&) {
            rejectedConcurrentOwner = true;
        }
        require(rejectedConcurrentOwner,
        "A second owner must not attach while the region is in use");
    }

    const int missingDescriptor = shm_open(name.c_str(), O_RDWR, 0600);
    require(missingDescriptor == -1 && errno == ENOENT,
            "Destroying the region must unlink its shared-memory object");
}

void testGncStateEstimatePublication() {
    aocs::GncState state = aocs::initializeGnc();
    sharedMemoryMap memory{};
    sensorMessage& input = memory.sensorData;
    const double t_J2000_s = 100.0;
    input.GPS_time_s = t_J2000_s + aocs::ipc::kGpsSecondsAtJ2000;
    input.aocs_time_J2000_s = t_J2000_s;
    input.has_new_GPS_measurement = 1;
    input.GPS_position_ECEF_m = {7000000.0, 1000.0, -2000.0};
    input.GPS_velocity_ECEF_m_s = {1.0, 7500.0, -2.0};
    for (std::size_t axis = 0; axis < aocs::ipc::kStateDimension; ++axis) {
        input.GPS_covariance_row_major[
            axis * aocs::ipc::kStateDimension + axis] = 1.0;
    }
    input.GPS_covariance_row_major[1] = 0.25;
    input.GPS_covariance_row_major[aocs::ipc::kStateDimension] = 0.25;

    require(aocs::runGncStep(state, memory) == 0,
            "GNC step must publish its initial GPS-based estimate");
    require(memory.navigationEstimate.valid == 1,
            "GNC must mark its published estimate valid");
    requireNear(memory.navigationEstimate.time_J2000_s, t_J2000_s, 1e-12,
        "Published estimate must use seconds since J2000");

    const aocs::environment::JulianDate date =
        aocs::environment::JulianDate::fromJ2000Seconds(t_J2000_s);
    const Eigen::Matrix3d ecefToEci =
        aocs::environment::EarthRotation(date, date)
            .reduce()
            .gcrsToItrs.transpose();
    const Eigen::Vector3d position_ecef(
        input.GPS_position_ECEF_m[0],
        input.GPS_position_ECEF_m[1],
        input.GPS_position_ECEF_m[2]);
    const Eigen::Vector3d velocity_ecef(
        input.GPS_velocity_ECEF_m_s[0],
        input.GPS_velocity_ECEF_m_s[1],
        input.GPS_velocity_ECEF_m_s[2]);
    const Eigen::Vector3d omega_ecef(0.0, 0.0, 7.2921150e-5);
    const Eigen::Vector3d expectedPosition = ecefToEci * position_ecef;
    const Eigen::Vector3d expectedVelocity =
        ecefToEci * (velocity_ecef + omega_ecef.cross(position_ecef));

    for (std::size_t axis = 0; axis < 3; ++axis) {
        requireNear(
            memory.navigationEstimate.position_ECI_m[axis],
            expectedPosition(static_cast<Eigen::Index>(axis)),
            1e-9,
            "Published position must be converted from ECEF to ECI");
        requireNear(
            memory.navigationEstimate.velocity_ECI_m_s[axis],
            expectedVelocity(static_cast<Eigen::Index>(axis)),
            1e-12,
            "Published velocity must include Earth-rotation conversion");
    }
    requireNear(memory.navigationEstimate.covariance_row_major[1], 0.25, 0.0,
        "Published covariance must preserve row-major measurement data");

    input.has_new_GPS_measurement = 0;
    input.aocs_time_J2000_s = t_J2000_s + 0.1;
    require(aocs::runGncStep(state, memory) == 0,
            "GNC must propagate between GPS measurements");
    requireNear(memory.navigationEstimate.time_J2000_s, t_J2000_s + 0.1,
        1e-12,
        "GNC must propagate its estimate at the AOCS timestep");
    const Eigen::Vector3d propagatedPosition(
        memory.navigationEstimate.position_ECI_m[0],
        memory.navigationEstimate.position_ECI_m[1],
        memory.navigationEstimate.position_ECI_m[2]);
    require((propagatedPosition - expectedPosition).norm() > 1.0,
            "GNC must propagate position between GPS updates");
    const double largestCovarianceMagnitude =
        Eigen::Map<const Eigen::Matrix<double, 6, 6, Eigen::RowMajor>>(
            memory.navigationEstimate.covariance_row_major.data())
            .cwiseAbs()
            .maxCoeff();
    require(largestCovarianceMagnitude < 10.0,
            "Gravity-gradient covariance propagation must remain bounded for a 0.1 s step");
    require(
        memory.navigationEstimate.covariance_row_major[3 * 6 + 3] > 1.0,
        "Acceleration process noise must contribute to propagated velocity covariance");

    input.aocs_time_J2000_s = t_J2000_s + 0.2;
    input.GPS_time_s =
        input.aocs_time_J2000_s + aocs::ipc::kGpsSecondsAtJ2000;
    input.has_new_GPS_measurement = 1;
    const double covarianceBeforeGpsUpdate =
        Eigen::Map<const Eigen::Matrix<double, 6, 6, Eigen::RowMajor>>(
            memory.navigationEstimate.covariance_row_major.data())(0, 0);
    require(aocs::runGncStep(state, memory) == 0,
            "GNC must process a subsequent GPS update");
    const double covarianceAfterGpsUpdate =
        Eigen::Map<const Eigen::Matrix<double, 6, 6, Eigen::RowMajor>>(
            memory.navigationEstimate.covariance_row_major.data())(0, 0);
    require(covarianceAfterGpsUpdate > 0.5 &&
                covarianceAfterGpsUpdate < covarianceBeforeGpsUpdate,
            "GPS update must reduce the propagated covariance without replacing the prior covariance");
}

}  // namespace

int main() {
    try {
        testQuaternion();
        testIau2000B();
        testRK4Propagation();
        testSharedMemoryGncTransport();
        testSharedMemoryRegionRecoversStaleObject();
        testGncStateEstimatePublication();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
