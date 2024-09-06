#ifndef MATH_H
#define MATH_H
#include <Eigen/Dense>

namespace math {
    bool floatEquals(double a, double b);
    constexpr double rad_to_deg = 57.29578;
    constexpr double deg_to_rad = 0.01745;
    double cot(double x);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::Matrix3d rotate_x(double radians);
    Eigen::Matrix3d rotate_y(double radians);
    Eigen::Matrix3d rotate_z(double radians);
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Vector3d euler_zyx_from_rotation(Eigen::Matrix3d &r);
    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v);
}
#endif
