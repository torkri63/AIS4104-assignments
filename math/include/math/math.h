#ifndef MATH_H
#define MATH_H
#include <Eigen/Dense>

namespace math {
    bool floatEquals(double a, double b);
    constexpr double rad_to_deg = 57.29578;
    constexpr double deg_to_rad = 0.01745;
    double cot(double x);
    Eigen::Matrix4d create_transformation_matrix(const Eigen::Matrix3d &R, const Eigen::Vector3d &p);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::Matrix3d rotate_x(double radians);
    Eigen::Matrix3d rotate_y(double radians);
    Eigen::Matrix3d rotate_z(double radians);
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e);
    Eigen::Vector3d euler_zyx_from_rotation(Eigen::Matrix3d &r);
    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v);
    Eigen::VectorXd sum_of_wrenches();
    Eigen::Matrix3d exponential_to_rotation_matrix(const Eigen::Vector3d &w, double theta);
    std::pair<Eigen::Vector3d, double> rotation_matrix_to_exponential(const Eigen::Matrix3d &r);
    Eigen::Matrix4d exponential_to_transformation_matrix(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);
    Eigen::Matrix4d exponential_to_transformation_matrix(const Eigen::VectorXd &screw, double theta);
    std::pair<Eigen::Vector3d, double> transformation_matrix_to_exponential(const Eigen::Matrix4d &T);
    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);
    Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions);
    Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions);
    }

#endif
