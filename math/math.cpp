#include <iostream>

#include <Eigen/Dense>
#include <math/math.h>
// Constants and equals
bool math::floatEquals(double a, double b)
{
    return std::abs(a - b) < 1e-6;
}
constexpr double rad_to_deg = 57.29578;
constexpr double deg_to_rad = 0.01745;

// Create a twist vector from velocity vectors
Eigen::VectorXd math::twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v) {
    Eigen::VectorXd twist(6);
    twist << w(0), w(1), w(2), v(0), v(1), v(2);
    return twist;
}

// Create a rotation matrix from rotating θ degrees about the principal axis x.
Eigen::Matrix3d math::rotate_x(double radians)
{
    Eigen::Matrix3d matrix;
    // Rotate about the x-axis, R_x
    matrix <<
        1, 0, 0,
        0, std::cos(radians), -std::sin(radians),
        0, std::sin(radians), std::cos(radians);
    return matrix;
}

// Create a rotation matrix from rotating θ degrees about the principal axis y.
Eigen::Matrix3d math::rotate_y(double radians)
{
    Eigen::Matrix3d matrix;
    // Rotate about the y-axis, R_y
    matrix <<
        std::cos(radians), 0, std::sin(radians),
        0, 1, 0,
        -std::sin(radians), 0, std::cos(radians);
    return matrix;
}

// Create a rotation matrix from rotating θ degrees about the principal axis z.
Eigen::Matrix3d math::rotate_z(double radians)
{
    Eigen::Matrix3d matrix;
    // Rotate about the z-axis, R_z
    matrix <<
        std::cos(radians), -std::sin(radians), 0,
        std::sin(radians), std::cos(radians), 0,
        0, 0, 1;
    return matrix;
}

// Create a rotation matrix from Euler angles (ZYX).
Eigen::Matrix3d math::rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    // Euler angles (IN RADIANS)
    double alpha = e[0]; // Rotation around Z-axis
    double beta = e[1];  // Rotation around Y-axis
    double gamma = e[2]; // Rotation around X-axis

    // Rotate around the individual matrices
    Eigen::Matrix3d Rz = rotate_z(alpha);
    Eigen::Matrix3d Ry = rotate_y(beta);
    Eigen::Matrix3d Rx = rotate_x(gamma);

    // Combine the rotations to a ZYX Euler rotation
    Eigen::Matrix3d rotation_matrix = Rz * Ry * Rx;
    return rotation_matrix;
}

// Compute Euler angles ZYX from rotation matrix
Eigen::Vector3d math::euler_zyx_from_rotation(Eigen::Matrix3d &r) {
    double a, b, c;
    if(floatEquals(r(2,0),-1.0)) {
        b = EIGEN_PI /2.0;
        a = 0.0;
        c = std::atan2(r(0,1), r(1,1));
    }
    else if(floatEquals(r(2,0),1.0)) {
        b = -EIGEN_PI/2.0;
        a = 0.0;
        c = -std::atan2(r(0,1), r(1,1));
    }
    else {
        b = std::atan2(-r(2,0), std::sqrt(r(0,0)*r(0,0) + r(1,0)*r(1,0)));
        a = std::atan2(r(1,0), r(0,0));
        c = std::atan2(r(2,1), r(2,2));
    }
    return Eigen::Vector3d(a, b, c);
}
















