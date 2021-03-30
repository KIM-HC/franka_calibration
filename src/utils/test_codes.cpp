// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

// test code ---------------------------------------------------------------------------
// header file

// cpp file

// // --------------------------------------------------------------------------- test code

// // char test code ---------------------------------------------------------------------------
// // header file
// #include <string>
// const char test_char[] = "/home/kimhc/git/robot_calibration/data/fixed_calibration/";

// // cpp file
// std::string test_string(test_char);
// test_string += " <- testing char";
// std::string test_better_string(std::string("testing better char: ") + std::string(test_char));
// std::cout << test_string << std::endl;
// std::cout << test_better_string << std::endl;
// std::string test_string_2 = std::string(test_char) + "<- testing char";
// // --------------------------------------------------------------------------- char test code


// // cpp-yaml test code ---------------------------------------------------------------------------
// // yaml file
// calibrate_base : true
// name : panda_left
// base : [0.0, 0.3, 1.0]

// relations:
//   - [1.0, 2.0, 3.0]
//   - [4.0, 5.0, 6.0]

// tests:
//   test1:
//     tt : 2
//   test2:
//     tt : 3
//   test21:
//     tt : 3

// // header file
// #include <string>
// #include <Eigen/Dense>
// #include "yaml-cpp/yaml.h"
// YAML::Node yaml_reader;

// // cpp file
// // check boolean
// if (yaml_reader["calibrate_base"]) {
//     std::cout << "succeeded reading boolean" << std::endl;
//     std::cout << "tester[\"calibrate_base\"]: " << yaml_reader["calibrate_base"] << std::endl;
// }

// // check eigen vector
// std::cout << "tester[\"base\"]: " << yaml_reader["base"] << std::endl;
// Eigen::Vector3d ev_(yaml_reader["base"].as<std::vector<double>>().data());
// std::cout << "eigen_vector base: " << ev_.transpose() << std::endl;

// // check string
// std::string yaml_name = yaml_reader["name"].as<std::string>();
// std::cout << "tester[\"name\"]: " << yaml_name << std::endl;

// // check others
// std::cout << "length relations: " << yaml_reader["relations"].size() << std::endl;
// std::cout << "tester[\"relations\"]: " << yaml_reader["relations"] << std::endl;
// std::cout << "tester[\"relations\"][0]: " << yaml_reader["relations"][0] << std::endl;
// std::cout << "tester[\"relations\"][0][1]: " << yaml_reader["relations"][0][1] << std::endl;

// std::cout << "length tests: " << yaml_reader["tests"].size() << std::endl;
// // --------------------------------------------------------------------------- cpp-yaml test code
