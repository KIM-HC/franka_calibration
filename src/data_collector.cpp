// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "robot_model/franka_panda_model.h"

#define DEBUG_MODE

std::string current_workspace = "/home/kimhc/git/robot_calibration/data/";

Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");

const int N_ROBOT = 2;  // number of robots
const int N_J = 7;    // number of joints in one robot
const int N_TOT = N_ROBOT * N_J;    // total joints
const int N_TOT_w_DATA = N_TOT + 1;  // total joints + data_number
const double DEL_MOVE = 0.001;

const int NUM_Q = N_TOT;

void forClosedChainCalibration(int data_number) {
  std::string data_input = current_workspace + "closed_chain_calibration/input_data/q_same_dir.txt";
  if (data_number == 2) {
    // 0.002
    data_input = current_workspace + "closed_chain_calibration/input_data/q_diff_dir.txt";
  } else if (data_number == 3) {
    // 0.004
    data_input = current_workspace + "closed_chain_calibration/input_data/q_same_dir.txt";
  }
#ifdef DEBUG_MODE
  std::cout << "reading from -- " << data_input << std::endl;
#endif
  std::string data_output = current_workspace + "closed_chain_calibration/input_data/input_data_" + std::to_string(data_number) + ".txt";
  int input_counter_ = 0;
  int output_counter_ = 0;
  Eigen::Matrix<double, N_TOT, 1> q_1;
  q_1.setZero();
  std::vector<Eigen::Matrix<double, NUM_Q, 1>> saver_;
  std::ifstream rf;

  rf.open(data_input);
  while (!rf.eof()) {
    input_counter_++;
    Eigen::Matrix<double, N_TOT, 1> new_q_1;
    for (int j=0; j < N_TOT; j++) {rf >> new_q_1(j);}
    if (N_ROBOT == 2) {
      Eigen::Matrix<double, N_J, 1> dump_;
      for (int j=0; j < N_J; j++) {rf >> dump_(j);}
    }

    if ((q_1-new_q_1).norm() > DEL_MOVE) {
      output_counter_++;
      Eigen::Matrix<double, NUM_Q, 1> saving_q;
      saving_q.head(N_TOT) = new_q_1;
      if (NUM_Q == N_TOT_w_DATA) {
        saving_q(N_TOT) = data_number;
        if (data_number == 3) {saving_q(N_TOT) = 1;}
      }
      saver_.push_back(saving_q);
      q_1 = new_q_1;
    }
  }
  rf.close();

  std::ofstream data_output_stream(data_output);
  for (int i=0; i < saver_.size(); i++) {
    if (i == saver_.size()-1)
      data_output_stream << saver_[i].transpose().format(tab_format);
    else
      data_output_stream << saver_[i].transpose().format(tab_format) << std::endl;
  }
  data_output_stream.close();
#ifdef DEBUG_MODE
  std::cout << "input size: " << input_counter_ << "\noutput size: " << output_counter_ << std::endl;
  std::cout << "saving at -- " << data_output << std::endl;
#endif
}

void for_fixed_calibration(int data_number) {
  std::string arm_name = "left";
  std::string data_input = current_workspace + "fixed_calibration/input_data/input_data_" + arm_name + "_" + std::to_string(data_number) + ".txt";
#ifdef DEBUG_MODE
  std::cout << "reading from -- " << data_input << std::endl;
#endif
  std::string data_output = current_workspace + "fixed_calibration/input_data/input_data_" + arm_name + "_" + std::to_string(data_number) + "_sorted.txt";
  int input_counter_ = 0;
  int output_counter_ = 0;
  Eigen::Matrix<double, 7, 1> q_1;
  q_1.setZero();
  std::vector<Eigen::Matrix<double, 7, 1>> saver_;
  std::ifstream rf;

  rf.open(data_input);
  while (!rf.eof()) {
    input_counter_++;
    Eigen::Matrix<double, 7, 1> new_q_1;
    for (int j=0; j < 7; j++) {rf >> new_q_1(j);}

    if ((q_1-new_q_1).norm() > DEL_MOVE) {
      output_counter_++;
      Eigen::Matrix<double, 7, 1> saving_q = new_q_1;
      saver_.push_back(saving_q);
      q_1 = new_q_1;
    }
  }
  rf.close();

  std::ofstream data_output_stream(data_output);
  for (int i=0; i < saver_.size(); i++) {
    if (i == saver_.size()-1)
      data_output_stream << saver_[i].transpose().format(tab_format);
    else
      data_output_stream << saver_[i].transpose().format(tab_format) << std::endl;
  }
  data_output_stream.close();
#ifdef DEBUG_MODE
  std::cout << "input size: " << input_counter_ << "\noutput size: " << output_counter_ << std::endl;
  std::cout << "saving at -- " << data_output << std::endl;
#endif
}

int main(int argc, char**argv) {
  int data_number = 1;
  if (argc >= 2) {
    data_number = strtod(argv[1], NULL);
  }
  forClosedChainCalibration(data_number);
  // for_fixed_calibration(data_number);
  return 0;
}
