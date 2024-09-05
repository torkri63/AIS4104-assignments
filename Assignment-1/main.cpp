#include <iostream>

#include <Eigen/Dense>
// *******************************************************************************
/* TASK 1: Skew symmetric matrix
   *******************************************************************************

    a) Create a function that calculates the skew-symmetric matrix representation of a vector v € R3.
    The function must be named skew_symmetric with return type Eigen::Matrix3d, and with one parameter
    of type Eigen::Vector3d named v.
    b) Verify the function using the function "skew_symmetric_test".
 */

// Function to calculate the skew-symmetric matrix of a vector v
Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& v) {
    Eigen::Matrix3d skewMat;
    skewMat <<    0, -v.z(),  v.y(),
               v.z(),     0, -v.x(),
              -v.y(),  v.x(),     0;
    return skewMat;
}

// Function to test the skew_symmetric function
void skew_symmetric_test() {
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});

    // Output the skew-symmetric matrix
    std::cout << "TASK 1:\n";
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;

    // Output the transposition of the skew-symmetric matrix
    std::cout << "Skew-symmetric matrix transposition (negated): " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

// *******************************************************************************
// TASK 2: Rotation matrices
// *******************************************************************************

//Create a rotation matrix from reference frame axes.
Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                const Eigen::Vector3d &y,
                                                const Eigen::Vector3d &z)
{
    Eigen::Matrix3d matrix;
    matrix.col(0) = x;
    matrix.col(1) = y;
    matrix.col(2) = z;
    return matrix;
}

// Create a rotation matrix from rotating θ degrees about the principal axis x.
Eigen::Matrix3d rotate_x(double degrees)
{
    Eigen::Matrix3d matrix;
    // Convert from degrees to Radians
    double radians = degrees * M_PI / 180;
    // Rotate about the x-axis, R_x
    matrix <<
        1, 0, 0,
        0, std::cos(radians), -std::sin(radians),
        0, std::sin(radians), std::cos(radians);
    return matrix;
}

// Create a rotation matrix from rotating θ degrees about the principal axis y.
Eigen::Matrix3d rotate_y(double degrees)
{
    Eigen::Matrix3d matrix;
    // Convert from degrees to Radians
    double radians = degrees * M_PI / 180;
    // Rotate about the y-axis, R_y
    matrix <<
        std::cos(radians), 0, std::sin(radians),
        0, 1, 0,
        -std::sin(radians), 0, std::cos(radians);
    return matrix;
}

// Create a rotation matrix from rotating θ degrees about the principal axis z.
Eigen::Matrix3d rotate_z(double degrees)
{
    Eigen::Matrix3d matrix;
    // Convert from degrees to Radians
    double radians = degrees * M_PI / 180;
    // Rotate about the z-axis, R_z
    matrix <<
        std::cos(radians), -std::sin(radians), 0,
        std::sin(radians), std::cos(radians), 0,
        0, 0, 1;
    return matrix;
}

// Create a rotation matrix from rotating θ degrees about an arbitrary axis.
Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees) {
    // Ensure the axis is a unit vector
    Eigen::Vector3d unit_axis = axis.normalized();

    // Convert degrees to radians
    double radians = degrees * M_PI / 180.0;

    // Define variables for Rodrigues rotation formula.
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(); // Identity matrix
    Eigen::Matrix3d u_skew = skew_symmetric(unit_axis); // Skew-symmetric matrix of the axis
    Eigen::Matrix3d u_skew_squared = u_skew * u_skew; // Square of the skew-symmetric matrix

    /* Input into Rodrigues rotation formula
     Rot(ω,θ)= I + sin θ[ω] + (1 - cos θ)[ω]^2*/
    Eigen::Matrix3d rotation_matrix = I + std::sin(radians) *u_skew + (1 - std::cos(radians)) * u_skew_squared;
    return rotation_matrix;
}

// Create a rotation matrix from Euler angles (ZYX).
Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    // Euler angles (in degrees)
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

// Verify the Functions for rotation matrices
void rotation_matrix_test()
{
    Eigen::Matrix3d rot =
    rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa =
    rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa =
    rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
                                    Eigen::Vector3d{-0.5, -0.5, 0.707107},
                                    Eigen::Vector3d{0.707107, -0.707107, 0.0});
    std::cout << "TASK 2:\n";
    std::cout << "Rotation matrix from Euler: " << std::endl;
    std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl;
    std::cout << "Rotation matrix from frame axes: " << std::endl;
    std::cout << rot_fa << std::endl << std::endl;
}


// *******************************************************************************
/* TASK 3: Transformation matrices
   *******************************************************************************

    Create a transformation matrix from a rotation matrix and translation vector.

    T = |R  p|
        |0  1|

    Verify the function for creating a transformation matrix
 */

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d matrix;
    // Set the top-left 3x3 block to the rotation matrix R
    matrix.block<3, 3>(0, 0) = r;

    // Set the top-right 3x1 block to the translation vector p
    matrix.block<3, 1>(0, 3) = p;

    // Set the bottom-left 1x3 block to zeroes
    matrix.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();

    // Set the bottom-right 1x1 block to 1
    matrix(3, 3) = 1.0;
    return matrix;
}

void transformation_matrix_test()
{
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "TASK 3:\n";
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl;
}

void transform_vector()
{
    // Define Euler angles which rotate frame alpha.
    Eigen::Vector3d euler_angles(60, 45, 0);

    // Define rotation matrix from Euler angles
    Eigen::Matrix3d R = rotation_matrix_from_euler_zyx(euler_angles);

    // Define p vector based on translation along z-axis
    Eigen::Vector3d p(0,0,10);

    // Define Homogenous Transformation Matrix based on R and p
    Eigen::Matrix4d T = transformation_matrix(R,p);

    // Define vector in Alpha frame
    Eigen::Vector3d V_alpha(2.5, 3.0, -10.0);

    // Transform Alpha vector to world frame
    Eigen::Vector4d V_alpha_H;
    V_alpha_H << V_alpha, 1;
    Eigen::Vector4d V_world_H = T * V_alpha_H;

    // Extract the transformed vector and print it
    Eigen::Vector3d vw = V_world_H.head<3>();
    std::cout << "Transformed vector in frame {w}:\n " << vw  << std::endl;
}






// Main function
int main()
{
    // Task 1
    skew_symmetric_test();
    //Task 2
    rotation_matrix_test();
    //
    transformation_matrix_test();
    transform_vector();
    return 0;
}















