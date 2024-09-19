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

// Create a rotation matrix from rotating about the principal axis x.
Eigen::Matrix3d math::rotate_x(double radians)
{
    Eigen::Matrix3d matrix;
    // Rotate about the x-axis, R_x
    matrix <<
        1.0, 0.0, 0.0,
        0.0, std::cos(radians), -std::sin(radians),
        0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

// Create a rotation matrix from rotating about the principal axis y.
Eigen::Matrix3d math::rotate_y(double radians)
{
    Eigen::Matrix3d matrix;
    // Rotate about the y-axis, R_y
    matrix <<
        std::cos(radians), 0.0, std::sin(radians),
        0.0, 1.0, 0.0,
        -std::sin(radians), 0.0, std::cos(radians);
    return matrix;
}

// Create a rotation matrix from rotating about the principal axis z.
Eigen::Matrix3d math::rotate_z(double radians)
{
    Eigen::Matrix3d matrix;
    // Rotate about the z-axis, R_z
    matrix <<
        std::cos(radians), -std::sin(radians), 0.0,
        std::sin(radians), std::cos(radians), 0.0,
        0.0, 0.0, 1.0;
    return matrix;
}

// Function to calculate the skew symmetric matrix of a vector v
Eigen::Matrix3d math::skew_symmetric(const Eigen::Vector3d &v)
{
    Eigen::Matrix3d skewMat;
    skewMat <<    0.0, -v.z(),  v.y(),
               v.z(),     0.0, -v.x(),
              -v.y(),  v.x(),     0.0;
    return skewMat;
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
    Eigen::Matrix3d R = Rz * Ry * Rx * Eigen::Matrix3d::Identity();
    return R;
}

// Transformation matrix from rotation and translation
Eigen::Matrix4d math::create_transformation_matrix(const Eigen::Matrix3d &R, const Eigen::Vector3d &p) {
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = p;
    return T;
}
// ASSIGNMENT 2

// Compute Euler angles ZYX from rotation matrix
Eigen::Vector3d math::euler_zyx_from_rotation(Eigen::Matrix3d &r) {
    double a, b, c;
    if(floatEquals(r(2, 0), -1.0))
    {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0, 1), r(1, 1));
    }
    else if(floatEquals(r(2, 0), 1.0))
    {
        b = -(EIGEN_PI / 2.0);
        a = 0.0;
        c = -std::atan2(r(0, 1), r(1, 1));
    }
    else
    {
        b = std::atan2(-r(2, 0), std::sqrt(r(0, 0) * r(0, 0) + r(1, 0) * r(1, 0)));
        a = std::atan2(r(1, 0), r(0, 0));
        c = std::atan2(r(2, 1), r(2, 2));
    }
    return Eigen::Vector3d{a, b, c};
}

// Create a twist vector from velocity vectors
Eigen::VectorXd math::twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v) {
    Eigen::VectorXd twist(6);
    twist << w(0), w(1), w(2), v(0), v(1), v(2);
    return twist;
}

// Create a screw axis from the axis of rotation, "s", point on the axis, "q" and pitch, "h".
Eigen::VectorXd math::screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h){
    Eigen::VectorXd screw_axis(6);
    Eigen::Vector3d v = -s.cross(q) + h * s;
    screw_axis << s, v;
    return screw_axis;
}

// Create the adjoint representation of a homogenous transformation matrix
Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf){
    // Define Adjoint matrix as 6x6
    Eigen::MatrixXd AdT(6, 6);

    // Fetch rotation matrix and translation vector from homogenous transformation matrix
    Eigen::Matrix3d R = tf.block<3, 3>(0, 0);
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);
    // Define zero matrix
    Eigen::Matrix3d zero = Eigen::Matrix3d::Zero();

    // Construct AdT matrix
    AdT << R,zero, skew_symmetric(p) * R, R;
    return AdT;
}

// Calculate the cotangent of x in RADIANS
double math::cot(double x){
    return std::cos(x) / std::sin(x); // Inverse of sin/cos = tan
}

// TASK 2

// Create rotation matrix from Euler angles YZX
Eigen::Matrix3d math::rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e) {
    // Euler angles (IN RADIANS)
    const double alpha = e[0]; // Rotation around Z-axis
    const double beta = e[1];  // Rotation around Y-axis
    const double gamma = e[2]; // Rotation around X-axis

    // Rotate around the individual matrices
    const Eigen::Matrix3d Rz = rotate_z(beta);
    const Eigen::Matrix3d Ry = rotate_y(alpha);
    const Eigen::Matrix3d Rx = rotate_x(gamma);

    // Combine the rotations to a ZYX Euler rotation
    Eigen::Matrix3d rotation_matrix = Ry * Rz * Rx;
    return rotation_matrix;
}

// Calculate sum of wrenches expressed in different reference frames
Eigen::VectorXd math::sum_of_wrenches() {
    // Define vectors and matrices
    Eigen::VectorXd Ff(6);
    Eigen::VectorXd Fh(6);
    Eigen::VectorXd Fa(6);
    Eigen::Matrix4d Thf;
    Eigen::Matrix4d Taf;

    // Assign values
    Fh << 0.0, 0.0, 0.0, 0.0, -5, 0.0;
    Fa << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    Thf <<
    1.0, 0.0, 0.0, -0.1,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0;
    Taf <<
    1.0, 0.0, 0.0, -0.25,
    0.0, 0.0, 1.0, 0.0,
    0.0, -1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0;
    // Calculate Ff
    Ff = adjoint_matrix(Thf).transpose() * Fh + adjoint_matrix(Taf).transpose() * Fa;
    return Ff;
}

// TASK 3

// Implement the matrix exponential for rotation matrices
Eigen::Matrix3d math::exponential_to_rotation_matrix(const Eigen::Vector3d &w, double theta){
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d w_hat = skew_symmetric(w);  // Skew-symmetric of vector omega

    // Rodrigues formula to compute rotation matrix
    Eigen::Matrix3d R = I + sin(theta) * w_hat + (1.0-cos(theta)) * w_hat*w_hat;
    return R;
}

// Implement matrix logarithm for rotation matrices
std::pair<Eigen::Vector3d, double> math::rotation_matrix_to_exponential(const Eigen::Matrix3d &R) {
    // Calculate the angle theta
    Eigen::AngleAxisd angle_axis(R); // Extracts the angle-axis representation from the rotation matrix
    Eigen::Vector3d w_hat = angle_axis.axis() * angle_axis.angle();
    double theta = angle_axis.angle();

    return std::make_pair(w_hat, theta);
}

Eigen::Matrix4d math::exponential_to_transformation_matrix(const Eigen::VectorXd &screw, double theta) {
    return exponential_to_transformation_matrix(
        {screw(0), screw(1), screw(2)},
        {screw(3), screw(4), screw(5)},
        theta);
}


// Implement the matrix exponential for homogeneous transformation matrices
Eigen::Matrix4d math::exponential_to_transformation_matrix(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta){
    // Fetch rotation matrix from earlier function
    const Eigen::Matrix3d w_hat = skew_symmetric(w);  // Skew-symmetric of vector omega
    const Eigen::Matrix3d R = exponential_to_rotation_matrix(w, theta);

    // Compute the translation part of T
    const Eigen::Vector3d p = (Eigen::Matrix3d::Identity()*theta + (1.0-cos(theta))*w_hat + (theta- sin(theta))*w_hat*w_hat) * v;

    // Construct the final transformation matrix
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = p;
    return T;
}

// Implement matrix logarithm for Transformation matrices
std::pair<Eigen::Vector3d, double> math::transformation_matrix_to_exponential(const Eigen::Matrix4d &T) {
    Eigen::Matrix3d R = T.block<3, 3>(0, 0);  // Extract rotation matrix
    Eigen::Vector3d p = T.block<3, 1>(0, 3);  // Extract translation vector

    // Compute the rotation angle theta
    double theta = std::acos((R.trace() - 1.0) / 2.0) * rad_to_deg;

    // Compute the skew-symmetric matrix omega_hat
    Eigen::Matrix3d omega_hat = (R - R.transpose()) / (2.0 * std::sin(theta*deg_to_rad));

    // Extract the rotation vector omega from the skew-symmetric omega_hat
    Eigen::Vector3d omega;
    omega << omega_hat(2, 1), omega_hat(0, 2), omega_hat(1, 0);

    // Compute the matrix V_inv
    Eigen::Matrix3d V_inv = Eigen::Matrix3d::Identity() - 0.5 * omega_hat + (1.0 / theta - 0.5 * std::tan((theta*deg_to_rad) / 2.0)) * (omega_hat * omega_hat);

    // Compute the translation part v
    Eigen::Vector3d v = V_inv * p;

    // Combine omega and v into a single vector S
    Eigen::VectorXd S(6);
    S << omega, v;

    return std::make_pair(S, theta);
}

// TASK 4

void math::print_pose(const std::string &label, const Eigen::Matrix4d &T) {
    // Extract the linear position
    Eigen::Vector3d p = T.block<3, 1>(0, 3);
    Eigen::Matrix3d R = T.block<3, 3>(0, 0);

    // Compute the Euler ZYX angles
    Eigen::Vector3d euler_angles = euler_zyx_from_rotation(R)*rad_to_deg;

    // Print the pose
    std::cout << label << ":\n";
    std::cout << "Linear position XYZ: [" << p.transpose() << "]\n";
    std::cout << "Euler ZYX angles: [" << euler_angles.transpose() << "]\n";
}

Eigen::Matrix4d math::planar_3r_fk_transform(const std::vector<double> &joint_positions) {
    constexpr double l1 = 10, l2 = 10, l3 = 10;

    // Extract joint positions
    double theta1 = joint_positions[0]*deg_to_rad;
    double theta2 = joint_positions[1]*deg_to_rad;
    double theta3 = joint_positions[2]*deg_to_rad;

    // Define homogenous transfer matrices
    Eigen::Matrix4d T01 = create_transformation_matrix(rotate_z(theta1), Eigen::Vector3d(0, 0, 0));
    Eigen::Matrix4d T12 = create_transformation_matrix(rotate_z(theta2), Eigen::Vector3d(l1, 0, 0));
    Eigen::Matrix4d T23 = create_transformation_matrix(rotate_z(theta3), Eigen::Vector3d(l2, 0, 0));
    Eigen::Matrix4d T34 = create_transformation_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d(l3, 0, 0));
    // Final transformation matrix
    Eigen::Matrix4d T04 = T01 * T12 * T23 * T34;
    return T04;
}

Eigen::Matrix4d math::planar_3r_fk_screw(const std::vector<double> &joint_positions) {
    constexpr double l1 = 10, l2 = 10, l3 = 10;

    // Define the screw axes for each joint
    const Eigen::VectorXd S1 = screw_axis({0.0,0.0,0.0}, {0.0,0.0,1.0}, 0.0);
    const Eigen::VectorXd S2 = screw_axis({l1,0.0,0.0}, {0.0,0.0,1.0}, 0.0);
    const Eigen::VectorXd S3 = screw_axis({l1 + l2,0.0,0.0}, {0.0,0.0,1.0}, 0.0);

    // Extract joint positions
    const double theta1 = joint_positions[0]*deg_to_rad;
    const double theta2 = joint_positions[1]*deg_to_rad;
    const double theta3 = joint_positions[2]*deg_to_rad;

    // Calculate the transformation matrices using the PoE formula
    const Eigen::Matrix4d T01 = exponential_to_transformation_matrix(S1, theta1);
    const Eigen::Matrix4d T12 = exponential_to_transformation_matrix(S2, theta2);
    const Eigen::Matrix4d T23 = exponential_to_transformation_matrix(S3, theta3);

    const Eigen::Matrix4d M = create_transformation_matrix(
        Eigen::Matrix3d::Identity(),
        {l1+l2+l3, 0.0, 0.0} );

    // Final transformation matrix
    Eigen::Matrix4d T04 = T01 * T12 * T23 * M;
    return T04;
}

// TASK 5

Eigen::Matrix4d math::ur3e_fk_screw(const std::vector<double> &joint_positions) {
    constexpr double h1 {0.15185}, l1 {-0.24355}, l2 {-0.2132}, h2 {0.08535},
    w1{-0.13105}, w2{-0.0921};

    // Define the screw axes for each joint of UR3e
    const Eigen::VectorXd S1 = screw_axis({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 0.0); // Joint 1
    const Eigen::VectorXd S2 = screw_axis({0.0, 0.0, h1}, {0.0, -1.0, 0.0}, 0.0); // Joint 2
    const Eigen::VectorXd S3 = screw_axis({l1, 0.0, h1}, {0.0, -1.0, 0.0}, 0.0); // Joint 3
    const Eigen::VectorXd S4 = screw_axis({l1 + l2, 0.0, h1}, {0.0, -1.0, 0.0}, 0.0); // Joint 4
    const Eigen::VectorXd S5 = screw_axis({l1 + l2, w1, 0}, {0.0, 0.0, -1.0}, 0.0); // Joint 5
    const Eigen::VectorXd S6 = screw_axis({l1 + l2, 0, h1 - h2}, {0.0, -1.0, 0.0}, 0.0); // Joint 6

    // Extract joint positions
    const double theta1 = joint_positions[0] * deg_to_rad;
    const double theta2 = joint_positions[1] * deg_to_rad;
    const double theta3 = joint_positions[2] * deg_to_rad;
    const double theta4 = joint_positions[3] * deg_to_rad;
    const double theta5 = joint_positions[4] * deg_to_rad;
    const double theta6 = joint_positions[5] * deg_to_rad;

    // Calculate the transformation matrices for each joint using the PoE formula
    const Eigen::Matrix4d T01 = exponential_to_transformation_matrix(S1, theta1);
    const Eigen::Matrix4d T12 = exponential_to_transformation_matrix(S2, theta2);
    const Eigen::Matrix4d T23 = exponential_to_transformation_matrix(S3, theta3);
    const Eigen::Matrix4d T34 = exponential_to_transformation_matrix(S4, theta4);
    const Eigen::Matrix4d T45 = exponential_to_transformation_matrix(S5, theta5);
    const Eigen::Matrix4d T56 = exponential_to_transformation_matrix(S6, theta6);

    // Construct M matrix for zero configuration
    Eigen::Matrix4d M; //ok
    M <<  1,  0, 0, l1+l2,
          0,  0,-1, w1+w2,
          0,  1, 0, h1-h2,
          0,  0, 0, 1;

    // Final transformation matrix from base to end-effector
    Eigen::Matrix4d T06 = T01 * T12 * T23 * T34 * T45 * T56 * M;
    return T06;
}

/*Eigen::Matrix4d math::ur3e_fk_transform(const std::vector<double> &joint_positions) {
    constexpr double h1 {0.15185}, l1 {-0.24355}, l2 {-0.2132}, h2 {0.08535},
    w1{-0.13105}, w2{-0.0921};

    // Extract joint positions
    const double theta1 = joint_positions[0] * deg_to_rad;
    const double theta2 = joint_positions[1] * deg_to_rad;
    const double theta3 = joint_positions[2] * deg_to_rad;
    const double theta4 = joint_positions[3] * deg_to_rad;
    const double theta5 = joint_positions[4] * deg_to_rad;
    const double theta6 = joint_positions[5] * deg_to_rad;

    // Define homogenous transfer matrices
    const Eigen::Matrix4d T01 = create_transformation_matrix(rotate_z(theta1), {0,0,h1});
    const Eigen::Matrix4d T12 = create_transformation_matrix(rotate_y(theta2), {0,w1,0});
    const Eigen::Matrix4d T23 = create_transformation_matrix(rotate_y(theta3), {l1,0,0});
    const Eigen::Matrix4d T34 = create_transformation_matrix(rotate_y(theta4), {l2,0,0});
    const Eigen::Matrix4d T45 = create_transformation_matrix(rotate_z(theta5), {0,0,-h2});
    const Eigen::Matrix4d T56 = create_transformation_matrix(rotate_y(theta6), {0,w2,0});

    Eigen::Matrix4d T06 = T01 * T12 * T23 * T34 * T45 * T56;
    return T06;
}*/

Eigen::Matrix4d math::ur3e_fk_transform(const std::vector<double> &joint_positions) {
    constexpr double h1 {0.15185}, l1 {-0.24355}, l2 {-0.2132}, h2 {0.08535},
                     w1 {-0.13105}, w2 {-0.0921};

    // Extract joint positions
    const double theta1 = joint_positions[0] * deg_to_rad;
    const double theta2 = joint_positions[1] * deg_to_rad;
    const double theta3 = joint_positions[2] * deg_to_rad;
    const double theta4 = joint_positions[3] * deg_to_rad;
    const double theta5 = joint_positions[4] * deg_to_rad;
    const double theta6 = joint_positions[5] * deg_to_rad;

    // Adjust rotation angles to match negative axes
    const Eigen::Matrix4d T01 = create_transformation_matrix(rotate_z(theta1),     {0, 0, h1});
    const Eigen::Matrix4d T12 = create_transformation_matrix(rotate_y(-theta2),    {0, w1, 0});
    const Eigen::Matrix4d T23 = create_transformation_matrix(rotate_y(-theta3),    {l1, 0, 0});
    const Eigen::Matrix4d T34 = create_transformation_matrix(rotate_y(-theta4),    {l2, 0, 0});
    const Eigen::Matrix4d T45 = create_transformation_matrix(rotate_z(-theta5),    {0, 0, -h2});
    const Eigen::Matrix4d T56 = create_transformation_matrix(rotate_y(-theta6),    {0, w2, 0});

    // Final transformation matrix from base to end-effector
    Eigen::Matrix4d T06 = T01 * T12 * T23 * T34 * T45 * T56;
    return T06;
}





