#pragma once

#include <Eigen/Dense>

Eigen::Matrix3d rot1(double angle_rad){
    Eigen::Matrix3d R;
    R << 1, 0, 0,
         0, cos(angle_rad), sin(angle_rad),
         0, -sin(angle_rad), cos(angle_rad);
    return R;
}

Eigen::Matrix3d rot2(double angle_rad){
    Eigen::Matrix3d R;
    R << cos(angle_rad), 0, -sin(angle_rad),
         0, 1, 0,
         sin(angle_rad), 0, cos(angle_rad);
    return R;
}

Eigen::Matrix3d rot3(double angle_rad){
    Eigen::Matrix3d R;
    R << cos(angle_rad), sin(angle_rad), 0,
         -sin(angle_rad), cos(angle_rad), 0,
         0, 0, 1;
    return R;
}