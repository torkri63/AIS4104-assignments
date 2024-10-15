find_package(Eigen3 CONFIG REQUIRED)

add_executable(Assignment-3 main.cpp)

target_link_libraries(Assignment-3 PUBLIC math PRIVATE Eigen3::Eigen)
