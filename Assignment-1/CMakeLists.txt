find_package(Eigen3 CONFIG REQUIRED)

add_executable(Assignment-1 main.cpp)

target_link_libraries(Assignment-1 PRIVATE Eigen3::Eigen)
