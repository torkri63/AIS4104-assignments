#include <iostream>
#include <math/math.h>

void ur3e_test_fk(){
std::cout << "Forward kinematics tests" << std::endl;
math::print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
math::print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
std::cout << std::endl;
math::print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0}) * math::deg_to_rad));
math::print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0}) * math::deg_to_rad));
std::cout << std::endl;
math::print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
math::print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
std::cout << std::endl;
math::print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
math::print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
}

void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0){
    auto [iterations, x_hat] = math::newton_raphson_root_find(f, x0);
    std::cout << "Newton Raphson root f, x0 = " << x0 << " -> iterations = " << iterations << " theta = " << x_hat << " f (theta) = " << f(x_hat) << std::endl;
}

void test_gradient_descent_root_find(const std::function<double(double)> &f, double x0){
    auto [iterations, x_hat] = math::gradient_descent_root_find(f, x0);
    std::cout << "Gradient Descent root f, x0 = " << x0 << " -> iterations = " << iterations << " theta = " << x_hat << " f (theta) = " << f(x_hat) << std::endl;
}

void test_root_find(){
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x)
    {
        return (x - 3.0) * (x - 3.0) - 1.0;
    };
    test_newton_raphson_root_find(f1, -20);
    test_gradient_descent_root_find(f1, -20);
}

void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions) {
    Eigen::Matrix4d tsb = math::ur3e_body_fk(joint_positions);
    auto [M, space_screws] = math::ur3e_space_chain();
    Eigen::MatrixXd jb = math::ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = math::ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = math::adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = math::adjoint_matrix(tsb.inverse());
    std::cout << "Jb: " << std::endl << jb << std::endl << "Ad_tbs*Js:" << std::endl << ad_tbs * js << std::endl << std::endl;
    std::cout << "Js: " << std::endl << js << std::endl << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb << std::endl << std::endl;
    std::cout << "d Jb: " << std::endl << jb - ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js - ad_tsb * jb << std::endl << std::endl;
}

void ur3e_test_jacobian(){
    std::cout << "Jacobian matrix tests" << std::endl;
    ur3e_test_jacobian(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad);
    ur3e_test_jacobian(math::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) * math::deg_to_rad);
}

void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0){
    std::cout << "Test from pose" << std::endl;
    Eigen::Matrix4d t_sd = math::create_transformation_matrix(math::rotation_matrix_from_euler_zyx(zyx), pos);
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = math::ur3e_body_fk(j_ik);
    math::print_pose(t_ik, " IK pose");
    math::print_pose(t_sd, "Desired pose");
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg << std::endl << std::endl;
}

void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0){
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = math::ur3e_space_fk(joint_positions);
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = math::ur3e_body_fk(j_ik);
    math::print_pose(t_ik, " IK pose");
    math::print_pose(t_sd, "Desired pose");
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg << std::endl << std::endl;
}

void ur3e_ik_test(){
    Eigen::VectorXd j_t0 = math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad;
    Eigen::VectorXd j_t1 = math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -89.0, 0.0, 0.0, 0.0}) * math::deg_to_rad;
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} * math::deg_to_rad, j_t0);
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} * math::deg_to_rad, j_t1);
    Eigen::VectorXd j_t2 = math::std_vector_to_eigen(std::vector<double>{50.0, -30.0, 20, 0.0, -30.0, 50.0}) * math::deg_to_rad;
    Eigen::VectorXd j_d1 = math::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) * math::deg_to_rad;
    ur3e_ik_test_configuration(j_d1, j_t0);
    ur3e_ik_test_configuration(j_d1, j_t2);
}


int main(){/*
    std::cout << "TASK 1" << std::endl;
    ur3e_test_fk();
    std::cout << "TASK 2" << std::endl;
    test_root_find();*/
    std::cout << "TASK 3" << std::endl;
    // ur3e_test_jacobian();
    std::cout << "TASK 4" << std::endl;
    ur3e_ik_test();
    return 0;
}





