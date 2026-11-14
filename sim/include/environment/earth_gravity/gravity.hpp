#pragma once

#include "environment/earth_gravity/egm2008_coefficients.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

enum class GravityModel {
    E_PointMass,
    E_EGM2008
};

class Gravity {
public:
    explicit Gravity(const EGM2008Coefficients* egm2008 = nullptr)
        : egm2008_(egm2008) {
        if (egm2008_ && (egm2008_->modelName != "EGM2008" ||
            egm2008_->norm != "fully_normalized" ||
            egm2008_->tideSystem != "tide_free")) {
            throw std::invalid_argument(
                "Gravity requires fully normalized, tide-free EGM2008 coefficients");
        }
    }

    Eigen::Vector3d computeGravity(
        const Eigen::Vector3d& position,
        double mass_kg,
        GravityModel gravity_model,
        std::uint32_t maxDegree,
        std::uint32_t maxOrder) const {
        const double radius = position.norm();
        if (!(radius > 0.0) || !(mass_kg > 0.0)) {
            throw std::invalid_argument("Position radius and spacecraft mass must be positive");
        }

        if (gravity_model == GravityModel::E_PointMass) {
            return -pointMassGM_ / (radius * radius * radius) * position * mass_kg;
        }
        if (gravity_model != GravityModel::E_EGM2008) {
            throw std::invalid_argument("Unsupported gravity model");
        }
        if (!egm2008_) {
            throw std::runtime_error("EGM2008 coefficients were not loaded");
        }
        return computeEgm2008(position, radius, maxDegree, maxOrder) * mass_kg;
    }

private:
    static double degreeRecurrenceA(std::uint32_t degree, std::uint32_t order) {
        const double n = degree;
        const double m = order;
        return std::sqrt((2.0 * n - 1.0) * (2.0 * n + 1.0) /
                         ((n - m) * (n + m)));
    }

    static double degreeRecurrenceB(std::uint32_t degree, std::uint32_t order) {
        const double n = degree;
        const double m = order;
        return std::sqrt((2.0 * n + 1.0) * (n - m - 1.0) * (n + m - 1.0) /
                         ((2.0 * n - 3.0) * (n - m) * (n + m)));
    }

    Eigen::Vector3d computeEgm2008(
        const Eigen::Vector3d& position,
        double radius,
        std::uint32_t maxDegree,
        std::uint32_t maxOrder) const {
        if (maxDegree > egm2008_->maxDegree || maxOrder > maxDegree) {
            throw std::invalid_argument(
                "EGM2008 order must not exceed degree, and degree must not exceed the dataset maximum");
        }

        const double sinLatitude = position.z() / radius;
        const double cosLatitude = std::hypot(position.x(), position.y()) / radius;
        const double longitude = std::atan2(position.y(), position.x());
        const double ratio = egm2008_->referenceRadius / radius;
        const double cosLongitude = std::cos(longitude);
        const double sinLongitude = std::sin(longitude);

        std::vector<double> cosOrder(maxOrder + 1);
        std::vector<double> sinOrder(maxOrder + 1);
        cosOrder[0] = 1.0;
        sinOrder[0] = 0.0;
        for (std::uint32_t order = 1; order <= maxOrder; ++order) {
            cosOrder[order] = cosOrder[order - 1] * cosLongitude -
                              sinOrder[order - 1] * sinLongitude;
            sinOrder[order] = sinOrder[order - 1] * cosLongitude +
                              cosOrder[order - 1] * sinLongitude;
        }

        std::vector<double> previous(maxOrder + 1, 0.0);
        std::vector<double> previousDerivative(maxOrder + 1, 0.0);
        std::vector<double> previous2(maxOrder + 1, 0.0);
        std::vector<double> previous2Derivative(maxOrder + 1, 0.0);
        std::vector<double> current(maxOrder + 1, 0.0);
        std::vector<double> currentDerivative(maxOrder + 1, 0.0);

        double radialSum = 0.0;
        double northSum = 0.0;
        double eastSum = 0.0;
        double degreeScale = 1.0;
        const double poleThreshold = 1.0e-12;
        const bool atPole = cosLatitude < poleThreshold;
        double pOverCosPrevious2 = 0.0;
        double pOverCosPrevious = std::sqrt(3.0);
        double pOverCosCurrent = 0.0;

        for (std::uint32_t degree = 0; degree <= maxDegree; ++degree) {
            if (degree == 0) {
                current[0] = 1.0;
                currentDerivative[0] = 0.0;
            } else if (degree == 1) {
                current[0] = std::sqrt(3.0) * sinLatitude;
                currentDerivative[0] = std::sqrt(3.0) * cosLatitude;
                if (maxOrder >= 1) {
                    current[1] = std::sqrt(3.0) * cosLatitude;
                    currentDerivative[1] = -std::sqrt(3.0) * sinLatitude;
                }
            } else {
                const std::uint32_t lastOrder = std::min(degree, maxOrder);
                for (std::uint32_t order = 0; order <= lastOrder; ++order) {
                    if (order == degree) {
                        const double factor = std::sqrt(
                            (2.0 * degree + 1.0) / (2.0 * degree));
                        current[order] = factor * cosLatitude * previous[order - 1];
                        currentDerivative[order] = factor *
                            (-sinLatitude * previous[order - 1] +
                             cosLatitude * previousDerivative[order - 1]);
                    } else {
                        const double factor = degreeRecurrenceA(degree, order);
                        current[order] = factor * sinLatitude * previous[order];
                        currentDerivative[order] = factor *
                            (cosLatitude * previous[order] +
                             sinLatitude * previousDerivative[order]);
                        if (degree >= order + 2) {
                            const double secondFactor = degreeRecurrenceB(degree, order);
                            current[order] -= secondFactor * previous2[order];
                            currentDerivative[order] -=
                                secondFactor * previous2Derivative[order];
                        }
                    }
                }
            }

            if (maxOrder == 0 || degree == 0) {
                pOverCosCurrent = 0.0;
            } else if (degree == 1) {
                pOverCosCurrent = std::sqrt(3.0);
            } else {
                const double factor = degreeRecurrenceA(degree, 1);
                pOverCosCurrent = factor * sinLatitude * pOverCosPrevious;
                if (degree >= 3) {
                    pOverCosCurrent -= degreeRecurrenceB(degree, 1) * pOverCosPrevious2;
                }
            }

            if (degree != 1) {
                const std::uint32_t lastOrder = std::min(degree, maxOrder);
                for (std::uint32_t order = 0; order <= lastOrder; ++order) {
                    const auto& coefficient = egm2008_->coefficient(
                        static_cast<std::uint16_t>(degree),
                        static_cast<std::uint16_t>(order));
                    const double c = coefficient.C;
                    const double s = coefficient.S;
                    const double cosine = cosOrder[order];
                    const double sine = sinOrder[order];
                    const double harmonic = c * cosine + s * sine;

                    radialSum += (degree + 1.0) * degreeScale *
                                 current[order] * harmonic;
                    northSum += degreeScale * currentDerivative[order] * harmonic;

                    if (order != 0) {
                        const double pOverCos = atPole
                            ? (order == 1 ? pOverCosCurrent : 0.0)
                            : current[order] / cosLatitude;
                        eastSum += degreeScale * order * pOverCos *
                                   (s * cosine - c * sine);
                    }
                }
            }

            previous2.swap(previous);
            previous.swap(current);
            previous2Derivative.swap(previousDerivative);
            previousDerivative.swap(currentDerivative);
            pOverCosPrevious2 = pOverCosPrevious;
            pOverCosPrevious = pOverCosCurrent;
            degreeScale *= ratio;
        }

        const double accelerationScale =
            egm2008_->earthGravityConstant / (radius * radius);
        const double radialAcceleration = -accelerationScale * radialSum;
        const double northAcceleration = accelerationScale * northSum;
        const double eastAcceleration = accelerationScale * eastSum;

        const Eigen::Vector3d radialDirection(
            cosLatitude * cosLongitude,
            cosLatitude * sinLongitude,
            sinLatitude);
        const Eigen::Vector3d northDirection(
            -sinLatitude * cosLongitude,
            -sinLatitude * sinLongitude,
            cosLatitude);
        const Eigen::Vector3d eastDirection(-sinLongitude, cosLongitude, 0.0);
        return radialAcceleration * radialDirection +
               northAcceleration * northDirection +
               eastAcceleration * eastDirection;
    }

    static constexpr double pointMassGM_ = 3.986004415e14; //using field as provided in coefficient file - see data/earth_gravity/README.md
    const EGM2008Coefficients* egm2008_;
};
