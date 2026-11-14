#include "environment/earth_rotation/iau2000b.hpp"
#include <math/quaternion.hpp>
#include "sim_constants.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

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

}  // namespace

int main() {
    try {
        testQuaternion();
        testIau2000B();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
