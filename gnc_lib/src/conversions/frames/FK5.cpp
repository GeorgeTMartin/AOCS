#include "math/matrix.hpp"
#include "conversions/frames/FK5.hpp"
#include "flight_constants.hpp"
#include <cmath>

/**
 * @brief Generates rotation quaternion from ECI to ECEF frame, given the time in seconds since J2000 epoch, using the FK5 model, ignoring polat motion and with truncated nutation.
 *
 * @param t_J2000TT_s Time in seconds since J2000 epoch.
 * @return Eigen::Quaterniond The corresponding quaternion.
 * 
 * Source:
 * [1] Vallado, David A. Fundamentals of Astrodynamics and Applications. 4th Ed.
 */

Eigen::Quaterniond getQuaternionECItoECEF(double t_J2000TT_s){
    Eigen::Quaterniond qECItoECEF = Eigen::Quaterniond::Identity();
    // Generate z, theta, zeta [1] Sec 3.7 Eqn. 3-88
    double t_J2000_julianCenturies = t_J2000TT_s / (86400.0 * 36525.0);

    double t_J2000_UT1_s = t_J2000TT_s + flight_constants::kTTtoUT1Offset_s;
    double t_J2000_UT1_julianCenturies = t_J2000_UT1_s / (86400.0 * 36525.0);

    double GMST_sec = fmod((67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * t_J2000_UT1_julianCenturies +
                  0.093104 * t_J2000_UT1_julianCenturies * t_J2000_UT1_julianCenturies -
                  6.2e-6 * t_J2000_UT1_julianCenturies * t_J2000_UT1_julianCenturies * t_J2000_UT1_julianCenturies), 86400.0);

    double GMST = (GMST_sec / 240.0) * flight_constants::kDegreesToRadians;
    

    double zeta = (2306.2181 * t_J2000_julianCenturies +
                    0.30188 * t_J2000_julianCenturies * t_J2000_julianCenturies +
                    0.017998 * t_J2000_julianCenturies * t_J2000_julianCenturies * t_J2000_julianCenturies) * flight_constants::kArcsecondsToRadians;

    double Theta = (2004.3109 * t_J2000_julianCenturies -
                    0.42665 * t_J2000_julianCenturies * t_J2000_julianCenturies -
                    0.041833 * t_J2000_julianCenturies * t_J2000_julianCenturies * t_J2000_julianCenturies) * flight_constants::kArcsecondsToRadians;

    double z = (2306.2181 * t_J2000_julianCenturies +
                1.09468 * t_J2000_julianCenturies * t_J2000_julianCenturies +
                0.018203 * t_J2000_julianCenturies * t_J2000_julianCenturies * t_J2000_julianCenturies) * flight_constants::kArcsecondsToRadians;

    double epsilon1980 = (84381.448 - 46.8150 * t_J2000_julianCenturies -
                        0.00059 * t_J2000_julianCenturies * t_J2000_julianCenturies +
                        0.001813 * t_J2000_julianCenturies * t_J2000_julianCenturies * t_J2000_julianCenturies) * flight_constants::kArcsecondsToRadians;

    double Omega = (125.04452 - 1934.136261 * t_J2000_julianCenturies) * flight_constants::kDegreesToRadians;

    double L = (280.4665 + 36000.7698 * t_J2000_julianCenturies) * flight_constants::kDegreesToRadians;
    double Lp = (218.3165 + 481267.8813 * t_J2000_julianCenturies) * flight_constants::kDegreesToRadians;

    double deltaEps = (9.20 * cos(Omega) + 0.57 * cos(2*L)
                 + 0.10 * cos(2*Lp) - 0.09 * cos(2*Omega)) * flight_constants::kArcsecondsToRadians;

    double deltaPsi = (-17.20 * sin(Omega) - 1.32 * sin(2*L)
                   - 0.23 * sin(2*Lp)  + 0.21 * sin(2*Omega)) * flight_constants::kArcsecondsToRadians;

    double epsilonTrue = epsilon1980 + deltaEps;

    double GAST = GMST + deltaPsi * cos(epsilon1980);

    // Generate PNRW matrices [1] Sec 3.7 Algorithm 24
    Eigen::Matrix3d precession = rot3(z) * rot2(-Theta) * rot3(zeta);
    Eigen::Matrix3d nutation = rot1(-epsilonTrue) * rot3(-deltaPsi) * rot1(epsilon1980);
    Eigen::Matrix3d rotation = rot3(GAST);

    Eigen::Matrix3d DCM_ECItoECEF = rotation * nutation * precession;
    qECItoECEF = Eigen::Quaterniond(DCM_ECItoECEF);



    return qECItoECEF;
}

Eigen::Quaterniond getQuaternionECEFtoECI(double t_J2000TT_s){
    Eigen::Quaterniond qECEFtoECI = getQuaternionECItoECEF(t_J2000TT_s).conjugate();
    return qECEFtoECI;
}