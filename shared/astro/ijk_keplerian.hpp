#pragma once

#include <Eigen/Dense>
#include "sim_constants.hpp"

/**
 * @brief Converts IJK coordinates to Keplerian elements.
 * 
 * Sources:
 * [1] Vallado, David A. Fundamentals of Astrodynamics and Applications. 4th ed.
 */

struct orbitalElements {
    double sma;
    double e;
    double i;
    double omega;
    double argp;
    double nu;
};
class IJKKeplerian {
public:
    static orbitalElements convert_ijk_to_keplerian(const Eigen::Vector3d& position_eci_m,const Eigen::Vector3d& velocity_eci_m_s) {
        // Source [1] Section 2.5 Algorithm 9
        Eigen::Vector3d vec3_AngularMomentum = position_eci_m.cross(velocity_eci_m_s);
        Eigen::Vector3d vec3_NodalVector = Eigen::Vector3d::UnitZ().cross(vec3_AngularMomentum);
        Eigen::Vector3d vec3_EccentricityVector = (velocity_eci_m_s.cross(vec3_AngularMomentum) / SimConstants::kEarthMu) - position_eci_m.normalized();
        
        double sma = -SimConstants::kEarthMu / (2.0 * (0.5 * velocity_eci_m_s.squaredNorm() - SimConstants::kEarthMu / position_eci_m.norm()));
        double e = vec3_EccentricityVector.norm();
        double i = std::acos(vec3_AngularMomentum.normalized().z());
        double omega = std::atan2(vec3_NodalVector.y(), vec3_NodalVector.x());
        
        const double small = 1.0e-10;
        const double r_norm = position_eci_m.norm();
        const double h_norm = vec3_AngularMomentum.norm();
        const double n_norm = vec3_NodalVector.norm();

        // Argument of perigee (argp)
        double argp = 0.0;
        if (n_norm > small && e > small) {
            // Using dot/cross ratio with proper scaling cancellation
            double cosine = vec3_NodalVector.dot(vec3_EccentricityVector) / (n_norm * e);
            double sine = (vec3_NodalVector.cross(vec3_AngularMomentum)).dot(vec3_EccentricityVector) / (n_norm * h_norm * e);
            argp = std::atan2(sine, cosine);
            if (argp < 0.0) argp += 2.0 * M_PI;
        }

        // True anomaly (nu)
        double nu = 0.0;
        if (e > small) {
            double cosine = vec3_EccentricityVector.dot(position_eci_m) / (e * r_norm);
            double sine = (vec3_EccentricityVector.cross(position_eci_m)).dot(vec3_AngularMomentum) / (e * h_norm * r_norm);
            nu = std::atan2(sine, cosine);
            if (nu < 0.0) nu += 2.0 * M_PI;
        } else {
            // Circular orbit true anomaly fallback (using argument of latitude)
            if (n_norm > small) {
                double cosine = vec3_NodalVector.dot(position_eci_m) / (n_norm * r_norm);
                double sine = (vec3_NodalVector.cross(vec3_AngularMomentum)).dot(position_eci_m) / (n_norm * h_norm * r_norm);
                nu = std::atan2(sine, cosine);
                if (nu < 0.0) nu += 2.0 * M_PI;
            } else {
                // Circular equatorial fallback
                nu = std::atan2(position_eci_m.y(), position_eci_m.x());
                if (nu < 0.0) nu += 2.0 * M_PI;
            }
        }

        return orbitalElements{sma, e, i, omega, argp, nu};
    }
};