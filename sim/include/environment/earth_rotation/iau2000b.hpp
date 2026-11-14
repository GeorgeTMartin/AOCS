#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include "sim_constants.hpp"

#include <math/quaternion.hpp>
#include <math/matrix/rotation.hpp>

namespace aocs::environment {

enum class EarthRotationModel {
    IAU2000B
};

struct JulianDate {
    double part1;
    double part2;

    static JulianDate fromJ2000Seconds(double secondsSinceJ2000) {
        return {2451545.0, secondsSinceJ2000 / 86400.0};
    }
};

struct NutationAngles {
    double longitude;
    double obliquity;
};

struct EarthRotationReduction {
    Eigen::Matrix3d gcrsToItrs;
    Eigen::Matrix3d gcrsToTrueEquinox;
    sim::math::Quaternion<double> gcrsToItrsQuaternion;
    NutationAngles nutation;
    double greenwichMeanSiderealTime;
    double equationOfEquinoxes;
    double greenwichApparentSiderealTime;
};

namespace detail {

constexpr double kTwoPi = 2.0 * SimConstants::kPi;
constexpr double kArcsecondsToRadians = SimConstants::kPi / (180.0 * 3600.0);
constexpr double kDaysPerCentury = 36525.0;
constexpr double kJ2000 = 2451545.0;

struct NutationTerm {
    int l;
    int lp;
    int f;
    int d;
    int om;
    double longitude;
    double longitudeRate;
    double longitudeCosine;
    double obliquity;
    double obliquityRate;
    double obliquitySine;
};

// Coefficients are in 0.1 microarcseconds, as in the IAU 2000B series.
// The coefficient tables are adapted from ERFA's nut00b.c and eect00.c;
// see ERFA_LICENSE.txt.
constexpr std::array<NutationTerm, 77> kNutationTerms{{
    {0, 0, 0, 0, 1, -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0},
    {0, 0, 2, -2, 2, -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0},
    {0, 0, 2, 0, 2, -2276413.0, -234.0, 2796.0, 978459.0, -485.0, 1374.0},
    {0, 0, 0, 0, 2, 2074554.0, 207.0, -698.0, -897492.0, 470.0, -291.0},
    {0, 1, 0, 0, 0, 1475877.0, -3633.0, 11817.0, 73871.0, -184.0, -1924.0},
    {0, 1, 2, -2, 2, -516821.0, 1226.0, -524.0, 224386.0, -677.0, -174.0},
    {1, 0, 0, 0, 0, 711159.0, 73.0, -872.0, -6750.0, 0.0, 358.0},
    {0, 0, 2, 0, 1, -387298.0, -367.0, 380.0, 200728.0, 18.0, 318.0},
    {1, 0, 2, 0, 2, -301461.0, -36.0, 816.0, 129025.0, -63.0, 367.0},
    {0, -1, 2, -2, 2, 215829.0, -494.0, 111.0, -95929.0, 299.0, 132.0},
    {0, 0, 2, -2, 1, 128227.0, 137.0, 181.0, -68982.0, -9.0, 39.0},
    {-1, 0, 2, 0, 2, 123457.0, 11.0, 19.0, -53311.0, 32.0, -4.0},
    {-1, 0, 0, 2, 0, 156994.0, 10.0, -168.0, -1235.0, 0.0, 82.0},
    {1, 0, 0, 0, 1, 63110.0, 63.0, 27.0, -33228.0, 0.0, -9.0},
    {-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0, 31429.0, 0.0, -75.0},
    {-1, 0, 2, 2, 2, -59641.0, -11.0, 149.0, 25543.0, -11.0, 66.0},
    {1, 0, 2, 0, 1, -51613.0, -42.0, 129.0, 26366.0, 0.0, 78.0},
    {-2, 0, 2, 0, 1, 45893.0, 50.0, 31.0, -24236.0, -10.0, 20.0},
    {0, 0, 0, 2, 0, 63384.0, 11.0, -150.0, -1220.0, 0.0, 29.0},
    {0, 0, 2, 2, 2, -38571.0, -1.0, 158.0, 16452.0, -11.0, 68.0},
    {0, -2, 2, -2, 2, 32481.0, 0.0, 0.0, -13870.0, 0.0, 0.0},
    {-2, 0, 0, 2, 0, -47722.0, 0.0, -18.0, 477.0, 0.0, -25.0},
    {2, 0, 2, 0, 2, -31046.0, -1.0, 131.0, 13238.0, -11.0, 59.0},
    {1, 0, 2, -2, 2, 28593.0, 0.0, -1.0, -12338.0, 10.0, -3.0},
    {-1, 0, 2, 0, 1, 20441.0, 21.0, 10.0, -10758.0, 0.0, -3.0},
    {2, 0, 0, 0, 0, 29243.0, 0.0, -74.0, -609.0, 0.0, 13.0},
    {0, 0, 2, 0, 0, 25887.0, 0.0, -66.0, -550.0, 0.0, 11.0},
    {0, 1, 0, 0, 1, -14053.0, -25.0, 79.0, 8551.0, -2.0, -45.0},
    {-1, 0, 0, 2, 1, 15164.0, 10.0, 11.0, -8001.0, 0.0, -1.0},
    {0, 2, 2, -2, 2, -15794.0, 72.0, -16.0, 6850.0, -42.0, -5.0},
    {0, 0, -2, 2, 0, 21783.0, 0.0, 13.0, -167.0, 0.0, 13.0},
    {1, 0, 0, -2, 1, -12873.0, -10.0, -37.0, 6953.0, 0.0, -14.0},
    {0, -1, 0, 0, 1, -12654.0, 11.0, 63.0, 6415.0, 0.0, 26.0},
    {-1, 0, 2, 2, 1, -10204.0, 0.0, 25.0, 5222.0, 0.0, 15.0},
    {0, 2, 0, 0, 0, 16707.0, -85.0, -10.0, 168.0, -1.0, 10.0},
    {1, 0, 2, 2, 2, -7691.0, 0.0, 44.0, 3268.0, 0.0, 19.0},
    {-2, 0, 2, 0, 0, -11024.0, 0.0, -14.0, 104.0, 0.0, 2.0},
    {0, 1, 2, 0, 2, 7566.0, -21.0, -11.0, -3250.0, 0.0, -5.0},
    {0, 0, 2, 2, 1, -6637.0, -11.0, 25.0, 3353.0, 0.0, 14.0},
    {0, -1, 2, 0, 2, -7141.0, 21.0, 8.0, 3070.0, 0.0, 4.0},
    {0, 0, 0, 2, 1, -6302.0, -11.0, 2.0, 3272.0, 0.0, 4.0},
    {1, 0, 2, -2, 1, 5800.0, 10.0, 2.0, -3045.0, 0.0, -1.0},
    {2, 0, 2, -2, 2, 6443.0, 0.0, -7.0, -2768.0, 0.0, -4.0},
    {-2, 0, 0, 2, 1, -5774.0, -11.0, -15.0, 3041.0, 0.0, -5.0},
    {2, 0, 2, 0, 1, -5350.0, 0.0, 21.0, 2695.0, 0.0, 12.0},
    {0, -1, 2, -2, 1, -4752.0, -11.0, -3.0, 2719.0, 0.0, -3.0},
    {0, 0, 0, -2, 1, -4940.0, -11.0, -21.0, 2720.0, 0.0, -9.0},
    {-1, -1, 0, 2, 0, 7350.0, 0.0, -8.0, -51.0, 0.0, 4.0},
    {2, 0, 0, -2, 1, 4065.0, 0.0, 6.0, -2206.0, 0.0, 1.0},
    {1, 0, 0, 2, 0, 6579.0, 0.0, -24.0, -199.0, 0.0, 2.0},
    {0, 1, 2, -2, 1, 3579.0, 0.0, 5.0, -1900.0, 0.0, 1.0},
    {1, -1, 0, 0, 0, 4725.0, 0.0, -6.0, -41.0, 0.0, 3.0},
    {-2, 0, 2, 0, 2, -3075.0, 0.0, -2.0, 1313.0, 0.0, -1.0},
    {3, 0, 2, 0, 2, -2904.0, 0.0, 15.0, 1233.0, 0.0, 7.0},
    {0, -1, 0, 2, 0, 4348.0, 0.0, -10.0, -81.0, 0.0, 2.0},
    {1, -1, 2, 0, 2, -2878.0, 0.0, 8.0, 1232.0, 0.0, 4.0},
    {0, 0, 0, 1, 0, -4230.0, 0.0, 5.0, -20.0, 0.0, -2.0},
    {-1, -1, 2, 2, 2, -2819.0, 0.0, 7.0, 1207.0, 0.0, 3.0},
    {-1, 0, 2, 0, 0, -4056.0, 0.0, 5.0, 40.0, 0.0, -2.0},
    {0, -1, 2, 2, 2, -2647.0, 0.0, 11.0, 1129.0, 0.0, 5.0},
    {-2, 0, 0, 0, 1, -2294.0, 0.0, -10.0, 1266.0, 0.0, -4.0},
    {1, 1, 2, 0, 2, 2481.0, 0.0, -7.0, -1062.0, 0.0, -3.0},
    {2, 0, 0, 0, 1, 2179.0, 0.0, -2.0, -1129.0, 0.0, -2.0},
    {-1, 1, 0, 1, 0, 3276.0, 0.0, 1.0, -9.0, 0.0, 0.0},
    {1, 1, 0, 0, 0, -3389.0, 0.0, 5.0, 35.0, 0.0, -2.0},
    {1, 0, 2, 0, 0, 3339.0, 0.0, -13.0, -107.0, 0.0, 1.0},
    {-1, 0, 2, -2, 1, -1987.0, 0.0, -6.0, 1073.0, 0.0, -2.0},
    {1, 0, 0, 0, 2, -1981.0, 0.0, 0.0, 854.0, 0.0, 0.0},
    {-1, 0, 0, 1, 0, 4026.0, 0.0, -353.0, -553.0, 0.0, -139.0},
    {0, 0, 2, 1, 2, 1660.0, 0.0, -5.0, -710.0, 0.0, -2.0},
    {-1, 0, 2, 4, 2, -1521.0, 0.0, 9.0, 647.0, 0.0, 4.0},
    {-1, 1, 0, 1, 1, 1314.0, 0.0, 0.0, -700.0, 0.0, 0.0},
    {0, -2, 2, -2, 1, -1283.0, 0.0, 0.0, 672.0, 0.0, 0.0},
    {1, 0, 2, 2, 1, -1331.0, 0.0, 8.0, 663.0, 0.0, 4.0},
    {-2, 0, 2, 2, 2, 1383.0, 0.0, -2.0, -594.0, 0.0, -2.0},
    {-1, 0, 0, 0, 2, 1405.0, 0.0, 4.0, -610.0, 0.0, 2.0},
    {1, 1, 2, -2, 2, 1290.0, 0.0, 0.0, -556.0, 0.0, 0.0}
}};

struct EquationComplementaryTerm {
    std::array<int, 8> multipliers;
    double sine;
    double cosine;
    double sineRate;
    double cosineRate;
};

constexpr std::array<EquationComplementaryTerm, 34> kEquationComplementaryTerms{{
    {{0, 0, 0, 0, 1, 0, 0, 0}, 2640.96e-6, -0.39e-6, 0.0, 0.0},
    {{0, 0, 0, 0, 2, 0, 0, 0}, 63.52e-6, -0.02e-6, 0.0, 0.0},
    {{0, 0, 2, -2, 3, 0, 0, 0}, 11.75e-6, 0.01e-6, 0.0, 0.0},
    {{0, 0, 2, -2, 1, 0, 0, 0}, 11.21e-6, 0.01e-6, 0.0, 0.0},
    {{0, 0, 2, -2, 2, 0, 0, 0}, -4.55e-6, 0.0, 0.0, 0.0},
    {{0, 0, 2, 0, 3, 0, 0, 0}, 2.02e-6, 0.0, 0.0, 0.0},
    {{0, 0, 2, 0, 1, 0, 0, 0}, 1.98e-6, 0.0, 0.0, 0.0},
    {{0, 0, 0, 0, 3, 0, 0, 0}, -1.72e-6, 0.0, 0.0, 0.0},
    {{0, 1, 0, 0, 1, 0, 0, 0}, -1.41e-6, -0.01e-6, 0.0, 0.0},
    {{0, 1, 0, 0, -1, 0, 0, 0}, -1.26e-6, -0.01e-6, 0.0, 0.0},
    {{1, 0, 0, 0, -1, 0, 0, 0}, -0.63e-6, 0.0, 0.0, 0.0},
    {{1, 0, 0, 0, 1, 0, 0, 0}, -0.63e-6, 0.0, 0.0, 0.0},
    {{0, 1, 2, -2, 3, 0, 0, 0}, 0.46e-6, 0.0, 0.0, 0.0},
    {{0, 1, 2, -2, 1, 0, 0, 0}, 0.45e-6, 0.0, 0.0, 0.0},
    {{0, 0, 4, -4, 4, 0, 0, 0}, 0.36e-6, 0.0, 0.0, 0.0},
    {{0, 0, 1, -1, 1, -8, 12, 0}, -0.24e-6, -0.12e-6, 0.0, 0.0},
    {{0, 0, 2, 0, 0, 0, 0, 0}, 0.32e-6, 0.0, 0.0, 0.0},
    {{0, 0, 2, 0, 2, 0, 0, 0}, 0.28e-6, 0.0, 0.0, 0.0},
    {{1, 0, 2, 0, 3, 0, 0, 0}, 0.27e-6, 0.0, 0.0, 0.0},
    {{1, 0, 2, 0, 1, 0, 0, 0}, 0.26e-6, 0.0, 0.0, 0.0},
    {{0, 0, 2, -2, 0, 0, 0, 0}, -0.21e-6, 0.0, 0.0, 0.0},
    {{0, 1, -2, 2, -3, 0, 0, 0}, 0.19e-6, 0.0, 0.0, 0.0},
    {{0, 1, -2, 2, -1, 0, 0, 0}, 0.18e-6, 0.0, 0.0, 0.0},
    {{0, 0, 0, 0, 0, 8, -13, -1}, -0.10e-6, 0.05e-6, 0.0, 0.0},
    {{0, 0, 0, 2, 0, 0, 0, 0}, 0.15e-6, 0.0, 0.0, 0.0},
    {{2, 0, -2, 0, -1, 0, 0, 0}, -0.14e-6, 0.0, 0.0, 0.0},
    {{1, 0, 0, -2, 1, 0, 0, 0}, 0.14e-6, 0.0, 0.0, 0.0},
    {{0, 1, 2, -2, 2, 0, 0, 0}, -0.14e-6, 0.0, 0.0, 0.0},
    {{1, 0, 0, -2, -1, 0, 0, 0}, 0.14e-6, 0.0, 0.0, 0.0},
    {{0, 0, 4, -2, 4, 0, 0, 0}, 0.13e-6, 0.0, 0.0, 0.0},
    {{0, 0, 2, -2, 4, 0, 0, 0}, -0.11e-6, 0.0, 0.0, 0.0},
    {{1, 0, -2, 0, -3, 0, 0, 0}, 0.11e-6, 0.0, 0.0, 0.0},
    {{1, 0, -2, 0, -1, 0, 0, 0}, 0.11e-6, 0.0, 0.0, 0.0},
    {{0, 0, 0, 0, 1, 0, 0, 0}, 0.0, 0.0, -0.87e-6, 0.0}
}};

double normalizeAngle(double angle) {
    angle = std::fmod(angle, kTwoPi);
    return angle < 0.0 ? angle + kTwoPi : angle;
}

double centuriesSinceJ2000(const JulianDate& date) {
    return ((date.part1 - kJ2000) + date.part2) / kDaysPerCentury;
}

Eigen::Matrix3d rotationX(double angle) {
    return MatrixRotation::rotationX(-angle);
}

Eigen::Matrix3d rotationY(double angle) {
    return MatrixRotation::rotationY(-angle);
}

Eigen::Matrix3d rotationZ(double angle) {
    return MatrixRotation::rotationZ(-angle);
}

double delaunayArgument(double constantArcseconds, double rateArcseconds,
                        double t) {
    return std::fmod(constantArcseconds + rateArcseconds * t, 1296000.0) *
           kArcsecondsToRadians;
}

NutationAngles nutationIau2000B(const JulianDate& tt) {
    const double t = centuriesSinceJ2000(tt);
    const double l = delaunayArgument(485868.249036, 1717915923.2178, t);
    const double lp = delaunayArgument(1287104.79305, 129596581.0481, t);
    const double f = delaunayArgument(335779.526232, 1739527262.8478, t);
    const double d = delaunayArgument(1072260.70369, 1602961601.2090, t);
    const double omega = delaunayArgument(450160.398036, -6962890.5431, t);

    double longitude = 0.0;
    double obliquity = 0.0;
    for (auto term = kNutationTerms.rbegin(); term != kNutationTerms.rend(); ++term) {
        const double argument = term->l * l + term->lp * lp + term->f * f +
                                term->d * d + term->om * omega;
        longitude += (term->longitude + term->longitudeRate * t) * std::sin(argument) +
                     term->longitudeCosine * std::cos(argument);
        obliquity += (term->obliquity + term->obliquityRate * t) * std::cos(argument) +
                     term->obliquitySine * std::sin(argument);
    }

    constexpr double kTenthMicroarcsecondToRadians =
        kArcsecondsToRadians / 1.0e7;
    constexpr double kMilliarcsecondToRadians =
        kArcsecondsToRadians / 1000.0;
    return {
        longitude * kTenthMicroarcsecondToRadians - 0.135 * kMilliarcsecondToRadians,
        obliquity * kTenthMicroarcsecondToRadians + 0.388 * kMilliarcsecondToRadians
    };
}

double equationComplementaryTerms(const JulianDate& tt) {
    const double t = centuriesSinceJ2000(tt);
    const std::array<double, 8> arguments{{
        delaunayArgument(485868.249036 + 31.8792 * t * t +
                         0.051635 * t * t * t - 0.00024470 * t * t * t * t,
                         1717915923.2178, t),
        delaunayArgument(1287104.793048 - 0.5532 * t * t +
                         0.000136 * t * t * t - 0.00001149 * t * t * t * t,
                         129596581.0481, t),
        delaunayArgument(335779.526232 - 12.7512 * t * t -
                         0.001037 * t * t * t + 0.00000417 * t * t * t * t,
                         1739527262.8478, t),
        delaunayArgument(1072260.703692 - 6.3706 * t * t +
                         0.006593 * t * t * t - 0.00003169 * t * t * t * t,
                         1602961601.2090, t),
        delaunayArgument(450160.398036 + 7.4722 * t * t +
                         0.007702 * t * t * t - 0.00005939 * t * t * t * t,
                         -6962890.5431, t),
        normalizeAngle(3.176146697 + 1021.3285546211 * t),
        normalizeAngle(1.753470314 + 628.3075849991 * t),
        normalizeAngle((0.024381750 + 0.00000538691 * t) * t)
    }};

    double sum = 0.0;
    for (auto term = kEquationComplementaryTerms.rbegin();
         term != kEquationComplementaryTerms.rend(); ++term) {
        double argument = 0.0;
        for (std::size_t index = 0; index < arguments.size(); ++index) {
            argument += term->multipliers[index] * arguments[index];
        }
        sum += (term->sine + term->sineRate * t) * std::sin(argument) +
               (term->cosine + term->cosineRate * t) * std::cos(argument);
    }
    return sum * kArcsecondsToRadians;
}

double earthRotationAngle(const JulianDate& ut1) {
    double first = ut1.part1;
    double second = ut1.part2;
    if (first > second) {
        std::swap(first, second);
    }
    const double daysSinceJ2000 = first + (second - kJ2000);
    const double fractionalDay = std::fmod(first, 1.0) + std::fmod(second, 1.0);
    return normalizeAngle(kTwoPi *
        (fractionalDay + 0.7790572732640 + 0.00273781191135448 * daysSinceJ2000));
}

void validateDate(const JulianDate& date) {
    if (!std::isfinite(date.part1) || !std::isfinite(date.part2)) {
        throw std::invalid_argument("Julian date parts must be finite");
    }
}

}  // namespace detail

class EarthRotation {
public:
    EarthRotation(JulianDate tt, JulianDate ut1, double xpRadians = 0.0,
                  double ypRadians = 0.0,
                  EarthRotationModel model = EarthRotationModel::IAU2000B)
        : tt_(tt), ut1_(ut1), xpRadians_(xpRadians), ypRadians_(ypRadians) {
        detail::validateDate(tt_);
        detail::validateDate(ut1_);
        if (!std::isfinite(xpRadians_) || !std::isfinite(ypRadians_)) {
            throw std::invalid_argument("Polar motion coordinates must be finite");
        }
        if (model != EarthRotationModel::IAU2000B) {
            throw std::invalid_argument("Unsupported Earth rotation model");
        }
    }

    EarthRotationReduction reduce() const {
        using namespace detail;

        const double t = centuriesSinceJ2000(tt_);
        const NutationAngles nutation = nutationIau2000B(tt_);
        const double meanObliquity =
            (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t) * t) * t) *
                kArcsecondsToRadians -
            0.02524 * t * kArcsecondsToRadians;

        const double dpsiBias = -0.041775 * kArcsecondsToRadians;
        const double depsBias = -0.0068192 * kArcsecondsToRadians;
        const double rightAscensionBias = -0.0146 * kArcsecondsToRadians;
        const double j2000Obliquity = 84381.448 * kArcsecondsToRadians;

        const double psi = (5038.7784 + (-1.07259 - 0.001147 * t) * t) *
                               t * kArcsecondsToRadians -
                           0.29965 * t * kArcsecondsToRadians;
        const double omega = j2000Obliquity +
            (0.05127 - 0.007726 * t) * t * t * kArcsecondsToRadians -
            0.02524 * t * kArcsecondsToRadians;
        const double chi = (10.5526 + (-2.38064 - 0.001125 * t) * t) *
                           t * kArcsecondsToRadians;

        const Eigen::Matrix3d bias =
            rotationX(-depsBias) *
            rotationY(dpsiBias * std::sin(j2000Obliquity)) *
            rotationZ(rightAscensionBias);
        const Eigen::Matrix3d precession =
            rotationZ(chi) * rotationX(-omega) *
            rotationZ(-psi) * rotationX(j2000Obliquity);
        const Eigen::Matrix3d nutationMatrix =
            rotationX(-(meanObliquity + nutation.obliquity)) *
            rotationZ(-nutation.longitude) * rotationX(meanObliquity);
        const Eigen::Matrix3d biasPrecessionNutation =
            nutationMatrix * precession * bias;

        const double era = earthRotationAngle(ut1_);
        const double gmstCentury = centuriesSinceJ2000(tt_);
        const double gmstOffsetArcseconds =
            0.014506 + (4612.15739966 +
            (1.39667721 + (-0.00009344 + 0.00001882 * gmstCentury) *
             gmstCentury) * gmstCentury) * gmstCentury;
        const double gmst = normalizeAngle(
            era + gmstOffsetArcseconds * kArcsecondsToRadians);

        const double equationOfEquinoxes =
            nutation.longitude * std::cos(meanObliquity) +
            equationComplementaryTerms(tt_);
        const double gast = normalizeAngle(gmst + equationOfEquinoxes);

        const double sp = -47.0e-6 * t * kArcsecondsToRadians;
        const Eigen::Matrix3d polarMotion =
            rotationX(-ypRadians_) * rotationY(-xpRadians_) * rotationZ(sp);
        const Eigen::Matrix3d gcrsToItrs =
            polarMotion * rotationZ(gast) * biasPrecessionNutation;

        return {
            gcrsToItrs,
            biasPrecessionNutation,
            sim::math::Quaternion<double>::fromPassiveMatrix(gcrsToItrs),
            nutation,
            gmst,
            equationOfEquinoxes,
            gast
        };
    }

private:
    JulianDate tt_;
    JulianDate ut1_;
    double xpRadians_;
    double ypRadians_;
};

}  // namespace aocs::environment
