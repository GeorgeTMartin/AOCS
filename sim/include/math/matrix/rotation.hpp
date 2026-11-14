#pragma once
#include "sim_constants.hpp"
#include <Eigen/Dense>

class MatrixRotation{
  public:
        static Eigen::Matrix3d rotationX(double angle) {
            return Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitX()).toRotationMatrix();
        }
        static Eigen::Matrix3d rotationY(double angle) {
            return Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitY()).toRotationMatrix();
        }
        static Eigen::Matrix3d rotationZ(double angle) {
            return Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        }
};