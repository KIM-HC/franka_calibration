// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#ifndef CLOSED_CHAIN_CALIBRATION_H_
#define CLOSED_CHAIN_CALIBRATION_H_

#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <vector>
#include <mutex>
#include <thread>
#include <string>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <Eigen/SVD>

#include "robot_model/franka_panda_model.h"
#include "utils/calibration_math.h"

const int N_ARM = 2;  // number of arms
const int N_PARAM = N_ARM * (N_ARM - 1);
const int N_DH = 4;   // number of dh parameters to calibrate
const int N_J = 7;    // number of joints in one robot
const int N_CAL = N_J * N_DH * N_ARM;  // total variables to calibrate
const int N_JDH = N_J * N_DH;    // number of dh parameters in one robot
const int PANDA_LEFT  =   0;
const int PANDA_RIGHT =   1;
const int PANDA_TOP   =   2;
int num_data;

// each vector's information is each arm's information
// each col (a, d, theta, alpha)
std::vector<Eigen::Matrix<double, N_J, N_DH>> offset_matrix;
// dataset
std::vector<Eigen::Matrix<double, N_ARM, N_J>> theta_data;
std::vector<Eigen::Isometry3d> T_W0;
std::vector<FrankaPandaModel> fpm;
std::vector<std::vector<double>> distTrue;
std::vector<int> data_number;

std::string  ws_ = "/home/kimhc/git/kinematics_calibration/";
std::string  data_ws_ = "data/closed_chain_calibration/";
std::string data_iter;
std::string data_offset;

Eigen::Isometry3d T_W0_LEFT;
Eigen::Isometry3d T_W0_RIGHT;
Eigen::Isometry3d T_W0_TOP;
Eigen::VectorXd del_phi(N_CAL);
Eigen::VectorXd p_total;
Eigen::MatrixXd jacobian;
Eigen::Matrix3d z_rot180;
Eigen::Matrix3d x_rot180;
Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");


#endif  // CLOSED_CHAIN_CALIBRATION_H_
