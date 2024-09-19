#include <iostream>
#include <math/math.h>

int main() {
    // Define euler angles, transform to rotation matrix and back
    Eigen::Vector3d e = Eigen::Vector3d(60.0,45.0,30.0) * math::deg_to_rad;
    Eigen::Matrix3d R = math::rotation_matrix_from_euler_zyx(e);
    Eigen::Vector3d ea = math::euler_zyx_from_rotation(R);
    std::cout << "\nTASK 1: \nImplementing algorithm, which receives euler zyx vector, transforms into rotation matrix and back to euler zyx vector" << std::endl;
    std::cout << "\nOriginal vector: [ " <<e.transpose() * math::rad_to_deg << "]" << std::endl;
    std::cout << "Vector after algorithm: [" <<ea.transpose() * math::rad_to_deg << "]" << std::endl;

    // Given data task 2
    std::cout << "\nTASK 2: Wrenches\n " << std::endl;
    Eigen::Vector3d fw(-30.0, 0.0, 0.0);  // Force in world frame
    Eigen::Vector3d ts(0.0, 0.0, 2.0);    // Torque in sensor frame
    Eigen::Vector3d ews = Eigen::Vector3d(60.0,-60.0,0.0) * math::deg_to_rad; // Euler angles
    // Calculate force in sensor frame and torque in world frame from rotation matrix
    Eigen::Matrix3d R_ws = math::rotation_matrix_from_euler_yzx(ews); // Define rotation matrix from Euler angles YZX
    Eigen::Vector3d fs = R_ws.transpose() * fw; // Transform force from world- to sensor frame
    Eigen::Vector3d tw = R_ws * ts; // Transform torque from sensor frame to world frame

    // Print vectors
    std::cout << "Force in world frame, f_w: [" << fw.transpose() << "] \n" << std::endl;
    std::cout << "Torque in world frame, t_w: [" << tw.transpose() << "] \n" << std::endl;
    std::cout << "Force in sensor frame, f_s: [" << fs.transpose() << "] \n" << std::endl;
    std::cout << "Torque in sensor frame, t_s: [" << ts.transpose() << "] \n" << std::endl;

    // Fetch sum of wrenches
    Eigen::VectorXd Ff = math::sum_of_wrenches();
    std::cout << "Sum of Wrenches;\nFf: [" << Ff.transpose() << "]" <<std::endl;

    // TASK 4
    std::cout << "\nTASK 4: \n" << std::endl;
    std::vector<std::vector<double>> test_positions = {
        {0.0, 0.0, 0.0},
        {90.0, 0.0, 0.0},
        {0.0, 90.0, 0.0},
        {0.0, 0.0, 90.0},
        {10.0, -15.0, 2.75}
    };
    std::cout << "\nForward kinematics using transformation matrices\n" << std::endl;
    for (const std::vector<double> &positions : test_positions) {
        Eigen::Matrix4d T = math::planar_3r_fk_transform(positions);
        std::cout << "----------------------------------------" << std::endl;
        math::print_pose("Pose for joint positions: [" + std::to_string(positions[0]) + ", " + std::to_string(positions[1]) + ", " + std::to_string(positions[2]) + "]", T);
    }
    std::cout << "\nForward kinematics using product of exponentials\n" << std::endl;
    for (const std::vector<double> &positions : test_positions) {
        Eigen::Matrix4d T = math::planar_3r_fk_screw(positions);
        std::cout << "----------------------------------------" << std::endl;
        math::print_pose("Pose for joint positions: [" + std::to_string(positions[0]) + ", " + std::to_string(positions[1]) + ", " + std::to_string(positions[2]) + "]", T);
    }

    std::cout << "\nTASK 5: \n" << std::endl;
    std::vector<std::vector<double>> test_positions_ur3e = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, -90.0, 0.0, 0.0},
        {0.0, -180.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -90.0, 0.0, 0.0, 0.0, 0.0}
    };

    std::cout << "\nForward kinematics for UR3e using Product of Exponentials\n" << std::endl;
    for (const std::vector<double> &positions : test_positions_ur3e) {
        Eigen::Matrix4d T = math::ur3e_fk_screw(positions);
        std::cout << "----------------------------------------" << std::endl;
        math::print_pose("Pose for joint positions: [" + std::to_string(positions[0]) + ", " + std::to_string(positions[1]) + ", " + std::to_string(positions[2]) + ", " + std::to_string(positions[3]) + ", " + std::to_string(positions[4]) + ", " + std::to_string(positions[5]) + "]", T);
    }
    std::cout << "\nForward kinematics for UR3e using Transformation Matrices\n" << std::endl;
    for (const std::vector<double> &positions : test_positions_ur3e) {
        Eigen::Matrix4d T = math::ur3e_fk_transform(positions);
        std::cout << "----------------------------------------" << std::endl;
        math::print_pose("Pose for joint positions: [" + std::to_string(positions[0]) + ", " + std::to_string(positions[1]) + ", " + std::to_string(positions[2]) + ", " + std::to_string(positions[3]) + ", " + std::to_string(positions[4]) + ", " + std::to_string(positions[5]) + "]", T);
    }
    return 0;
}





