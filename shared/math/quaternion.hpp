#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace sim::math {

// Scalar-first Hamilton quaternion. Rotation methods interpret q as a
// passive coordinate transform: v' = conjugate(q) * (0,v) * q. Therefore
// M(q1*q2) = M(q2)*M(q1) for the passive matrices returned below.
template <typename Scalar>
struct Quaternion {
    static_assert(std::is_floating_point<Scalar>::value,
                  "Quaternion scalar type must be floating point");

    Scalar w{1};
    Scalar x{0};
    Scalar y{0};
    Scalar z{0};

    constexpr Quaternion() = default;
    constexpr Quaternion(Scalar scalar, Scalar xComponent, Scalar yComponent,
                        Scalar zComponent)
        : w(scalar), x(xComponent), y(yComponent), z(zComponent) {}

    static constexpr Quaternion identity() {
        return {};
    }

    Scalar squaredNorm() const {
        return w * w + x * x + y * y + z * z;
    }

    Scalar norm() const {
        return std::sqrt(squaredNorm());
    }

    Quaternion normalized() const {
        const Scalar magnitude = norm();
        if (!(magnitude > Scalar(0)) || !std::isfinite(magnitude)) {
            throw std::invalid_argument("Cannot normalize a zero or non-finite quaternion");
        }
        return {w / magnitude, x / magnitude, y / magnitude, z / magnitude};
    }

    constexpr Quaternion conjugate() const {
        return {w, -x, -y, -z};
    }

    Quaternion inverse() const {
        const Scalar magnitudeSquared = squaredNorm();
        if (!(magnitudeSquared > Scalar(0)) || !std::isfinite(magnitudeSquared)) {
            throw std::invalid_argument("Cannot invert a zero or non-finite quaternion");
        }
        const Quaternion conjugated = conjugate();
        return {
            conjugated.w / magnitudeSquared,
            conjugated.x / magnitudeSquared,
            conjugated.y / magnitudeSquared,
            conjugated.z / magnitudeSquared
        };
    }

    constexpr Quaternion operator*(const Quaternion& rhs) const {
        return {
            w * rhs.w - x * rhs.x - y * rhs.y - z * rhs.z,
            w * rhs.x + x * rhs.w + y * rhs.z - z * rhs.y,
            w * rhs.y - x * rhs.z + y * rhs.w + z * rhs.x,
            w * rhs.z + x * rhs.y - y * rhs.x + z * rhs.w
        };
    }

    static Quaternion fromPassiveAxisAngle(
        const Eigen::Matrix<Scalar, 3, 1>& axis, Scalar angleRadians) {
        const Scalar axisNorm = axis.norm();
        if (!(axisNorm > Scalar(0)) || !std::isfinite(axisNorm) ||
            !std::isfinite(angleRadians)) {
            throw std::invalid_argument("Axis and angle must be finite and the axis nonzero");
        }
        const Scalar halfAngle = angleRadians / Scalar(2);
        const Scalar sine = std::sin(halfAngle) / axisNorm;
        return {std::cos(halfAngle), axis.x() * sine, axis.y() * sine, axis.z() * sine};
    }

    static Quaternion fromPassiveMatrix(
        const Eigen::Matrix<Scalar, 3, 3>& passiveMatrix) {
        if (!passiveMatrix.allFinite()) {
            throw std::invalid_argument("Rotation matrix must be finite");
        }

        const Eigen::Matrix<Scalar, 3, 3> identityMatrix =
            Eigen::Matrix<Scalar, 3, 3>::Identity();
        const Scalar tolerance = Scalar(1e-6);
        if ((passiveMatrix.transpose() * passiveMatrix - identityMatrix).norm() > tolerance ||
            std::abs(passiveMatrix.determinant() - Scalar(1)) > tolerance) {
            throw std::invalid_argument("Matrix must be a proper orthogonal rotation");
        }

        const Eigen::Quaternion<Scalar> active(passiveMatrix.transpose());
        return Quaternion{active.w(), active.x(), active.y(), active.z()}.normalized();
    }

    Eigen::Matrix<Scalar, 3, 1> rotatePassive(
        const Eigen::Matrix<Scalar, 3, 1>& vector) const {
        const Quaternion unit = normalized();
        const Quaternion pureVector(Scalar(0), vector.x(), vector.y(), vector.z());
        const Quaternion rotated = unit.conjugate() * pureVector * unit;
        return {rotated.x, rotated.y, rotated.z};
    }

    Eigen::Matrix<Scalar, 3, 3> toPassiveMatrix() const {
        const Quaternion unit = normalized();
        const Scalar ww = unit.w * unit.w;
        const Scalar xx = unit.x * unit.x;
        const Scalar yy = unit.y * unit.y;
        const Scalar zz = unit.z * unit.z;
        Eigen::Matrix<Scalar, 3, 3> result;
        result << ww + xx - yy - zz, Scalar(2) * (unit.x * unit.y + unit.w * unit.z),
                  Scalar(2) * (unit.x * unit.z - unit.w * unit.y),
                  Scalar(2) * (unit.x * unit.y - unit.w * unit.z), ww - xx + yy - zz,
                  Scalar(2) * (unit.y * unit.z + unit.w * unit.x),
                  Scalar(2) * (unit.x * unit.z + unit.w * unit.y),
                  Scalar(2) * (unit.y * unit.z - unit.w * unit.x), ww - xx - yy + zz;
        return result;
    }
};

}  // namespace sim::math
