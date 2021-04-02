// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#pragma once
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <vector>  // std::vector
#include <mutex>
#include <thread>
#include <string>  // std::string
#include <utility>  // std::pair

#include "yaml-cpp/yaml.h"
#include "robot_model/franka_panda_model.h"
#include "utils/calibration_math.h"

#define DEBUG_MODE

using calib_math::computePandaTransform;

struct RobotRelation {
    std::string base_arm;
    std::string target_arm;
    Eigen::Vector3d relation;  // from base to target
};
struct RobotState {
    Eigen::Isometry3d base;  // robot base
    std::string name;  // robot name
    int qp;  // place in q inputs
    FrankaPandaModel fpm;
};

Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");
const int num_dh_ {4};   // number of dh parameters to calibrate
const int num_j_ = 7;    // number of joints in one robot
const int per_rel_ {3};  // function per relation

class ClosedCalibration {
 public:
    // initialization part
    explicit ClosedCalibration(const std::string &yn);

    // exectues calibration
    void executeCalibration();

 private:
    // constants
    const std::string arm_priority[3] = {"panda_left", "panda_right", "panda_top"};
    const std::string ws_ {"/home/kimhc/git/robot_calibration/data/closed_chain_calibration/"};
    const std::string yws_ {ws_ + "yaml/"};
    const std::string iws_ {ws_ + "input_data/"};
    const std::string rws_ {ws_ + "result/"};
    const std::string dws_ {ws_ + "debug/"};

    // variables
    std::vector<std::vector<RobotRelation>> rr_;  // vec1: type | vec1: relations
    std::vector<std::string> file_names_;  // vec: type | string: file name
    std::vector<std::vector <Eigen::VectorXd>> inputs;  // vec1: type | vec2: inputs | eigen: q value
    std::vector<int> input_size_stack_;
    std::map<std::string, RobotState> rs_;  // robots
    std::map<int, std::string> rnn_;  // maps robot number and name
    std::map<std::string, Eigen::Matrix<double, 7, 1>> qv_;  // robots  // error - num_j_

    std::ofstream wfi_;  // iteration information
    std::ifstream rf_;  // reader

    std::string midpoint_path;  // path to midpoint save
    std::string yaml_path_;  // path to yaml
    std::string method_;  // types of method - lm_method
    std::string save_name_;  
    double lambda_;  // damping parameter for LM method
    double lambda_factor_ {1.1};  //
    bool continue_calib_;  // if true, continue from previous work
    bool calibrate_base_;  // if true, calibrate base also
    int num_rel_type_;  // number of types
    int num_rel_;  // number of relation
    int num_arm_;  // number of arm
    int tot_j_;  // total number of joints
    int sn_;  // save number
    int num_input_ {0};  // number of inputs (q data)

    Eigen::MatrixXd jacobian_, dh_mat_;
    Eigen::VectorXd del_pos_, del_dh_;
    Eigen::Matrix<double, 7, 4> panda_dh;

    // reads yaml
    void readYaml();

    // reads data & makes structure
    void structData();

    // loads data when continue_calib_ is true
    void loadData();

    // updates fpm and others if needed
    void updateStatus();

    // computes jacobian using analytical method
    void computeJacobian();

    // computes jacobian using analytical method
    Eigen::VectorXd computeDeltaPos();
};
