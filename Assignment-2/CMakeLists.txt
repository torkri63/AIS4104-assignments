find_package(Eigen3 CONFIG REQUIRED)

add_executable(Assignment-2 main.cpp
        ../Assignment-3/main.cpp)

target_link_libraries(Assignment-2 PUBLIC math PRIVATE Eigen3::Eigen)
