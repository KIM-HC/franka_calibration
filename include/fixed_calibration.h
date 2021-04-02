// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>

#include "yaml-cpp/yaml.h"
#include "robot_model/franka_panda_model.h"
#include "utils/calibration_math.h"

using calib_math::transformDH;

// #define TEST_PRINT

const int N_DH = 4;  // number of dh parameters to calibrate
const int N_J = 7;  // number of joints to calibrate
const int N_CAL = N_DH*N_J;  // total variables to calibrate

double lambda = 1.0;
int sn;
const char ws_[] {"/home/kimhc/git/robot_calibration/data/fixed_calibration/"};
Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");
Eigen::Matrix<double, 7, 4> panda_dh;
std::ofstream iteration_info;
std::ofstream mid_point_save;

// data set
struct dataSetStruct {
    std::string method;  // lm_method, svd_method
    std::string arm_name;
    std::string yaml_path;
    int num_ref;
    int num_input;
    bool calibrate_base;

    std::vector<Eigen::Vector3d> relations;
    std::vector<std::string> file_names;
    std::vector <std::pair<Eigen::Matrix<double, N_J, 1>, Eigen::Matrix<double, 3, 1>>> inputs;

    Eigen::Matrix<double, N_CAL, 1> dh_vec;
    Eigen::Matrix<double, N_CAL, 1> del_dh;
    Eigen::Matrix<double, N_J, N_DH> dh_mat;
    Eigen::Isometry3d base;
    Eigen::MatrixXd jacobian;
    Eigen::VectorXd del_pos;
};

// reads yaml
void readYaml(dataSetStruct *ds);

// reads data & makes structure
void structData(dataSetStruct *ds);

// computes jacobian using analytical method
void computeJacobian(dataSetStruct *ds, FrankaPandaModel *fpm);

// computes jacobian using analytical method
void computeDeltaPos(dataSetStruct *ds, FrankaPandaModel *fpm);

// writes pos info for debugging
void writePosInfo(dataSetStruct *ds, FrankaPandaModel *fpm);
