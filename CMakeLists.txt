cmake_minimum_required(VERSION 2.8.3)
project(kinematics_calibration)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(RBDL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories(
 include
 ${EIGEN3_INCLUDE_DIRS}    
 ${RBDL_INCLUDE_DIR}
 ${YAML_CPP_INCLUDE_DIR}
)

add_library(calibration_math_lib
  src/utils/calibration_math.cpp
)

target_link_libraries(calibration_math_lib
  ${catkin_LIBRARIES}
  ${RBDL_LIBRARY}
)

add_executable(fixed_calibration
  src/fixed_calibration.cpp
  src/robot_model/franka_panda_model.cpp
)
target_link_libraries(fixed_calibration
  ${YAML_CPP_LIBRARIES}
  calibration_math_lib
  pthread
)

add_executable(closed_chain_calibration
  src/closed_chain_calibration.cpp
  src/robot_model/franka_panda_model.cpp
)
target_link_libraries(closed_chain_calibration
  ${YAML_CPP_LIBRARIES}
  calibration_math_lib
  pthread
)

add_executable(data_collector
  src/data_collector.cpp
  src/robot_model/franka_panda_model.cpp
)
target_link_libraries(data_collector
  calibration_math_lib
  pthread
)
