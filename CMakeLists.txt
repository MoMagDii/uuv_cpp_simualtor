cmake_minimum_required(VERSION 3.5)
project(simulation)

add_definitions(-std=c++17)
# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++17
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)
find_package(ruckig REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Matplot++ REQUIRED)
add_subdirectory(libs/lin_alg_tools)
find_package(Python3 COMPONENTS NumPy)

add_executable(${PROJECT_NAME} src/main.cpp src/optimal_trajectory.cpp src/plant.cpp src/lqr.cpp)

target_link_libraries(${PROJECT_NAME} ${EIGEN3_LIBRARIES} ${LAPACK_LIBRARIES} ruckig ipopt Matplot++::matplot   Python3::Python Python3::Module)
target_include_directories(${PROJECT_NAME}
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/}
  ${EIGEN3_INCLUDE_DIRS}
  libs/lin_alg_tools/include)

