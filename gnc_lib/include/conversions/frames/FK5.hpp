#pragma once

#include <Eigen/Dense>

Eigen::Quaterniond getQuaternionECItoECEF(double time_J2000_s);
Eigen::Quaterniond getQuaternionECEFtoECI(double time_J2000_s);