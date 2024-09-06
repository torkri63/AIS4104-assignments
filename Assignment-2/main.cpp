#include <iostream>
#include <math/math.h>

int main() {
    Eigen::Vector3d e = Eigen::Vector3d(60,45,30) * math::deg_to_rad;
    Eigen::Matrix3d r = math::rotation_matrix_from_euler_zyx(e);
    Eigen::Vector3d ea = math::euler_zyx_from_rotation(r);
    std::cout << "Original vector: [ " <<e.transpose() * math::rad_to_deg << "]" << std::endl;
    std::cout << "Vector after algorithm: [" <<ea.transpose() * math::rad_to_deg << "]" << std::endl;
    return 0;
}






