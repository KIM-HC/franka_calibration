// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#include "../include/closed_chain_calibration.h"

int main(int argc, char**argv) {
  std::string yaml_name_;
  if (argc == 2) {
    yaml_name_ = argv[1];
  } else {
    yaml_name_ = "closed_1";
  }
  std::cout << "yaml name: " << yaml_name_ << std::endl;
  ClosedCalibration cc(yaml_name_);
  cc.executeCalibration();
  return 0;
}

ClosedCalibration::ClosedCalibration(const std::string &yn) {
  yaml_path_ = yws_ + yn + ".yaml";
  readYaml();
  structData();
  // if (continue_calib_)
  //   loadData();
}

void ClosedCalibration::readYaml() {
  YAML::Node yn_ = YAML::LoadFile(yaml_path_);

  num_rel_ = yn_["relation_types"][0]["relations"].size();
  num_rel_type_ = yn_["relation_types"].size();
  num_arm_ = yn_["robots"].size();
  calibrate_base_ = yn_["calibrate_base"];
  continue_calib_ = yn_["continue_calib"];
  method_ = yn_["method"].as<std::string>();
  save_name_ = yn_["save_name"].as<std::string>();
  lambda_ = yn_["lambda"].as<double>();

  std::vector<std::string> nl_;  // name list
  for (int i=0; i < num_arm_; i++) {
    RobotState st_;
    st_.name = yn_["robots"][i]["name"].as<std::string>();
    nl_.push_back(st_.name);
    st_.base.setIdentity();
    st_.base.translation() = Eigen::Vector3d(yn_["robots"][i]["base"].as<std::vector<double>>().data());
    st_.base.linear() = Eigen::AngleAxisd(yn_["robots"][i]["rot"].as<double>()*M_PI/180.0, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    rs_.insert(std::make_pair(st_.name, st_));
    std::cout << st_.name << " base:\n" << st_.base.matrix() << std::endl;
  }
  int qps_ = 0;
  for (int i=0; i < 3; i++) {
    auto it = std::find(nl_.begin(), nl_.end(), arm_priority[i]);
    if (it != nl_.end()) {
      #ifdef DEBUG_MODE
      std::cout << arm_priority[i] << " will get number " << qps_ << std::endl;
      #endif
      rs_[arm_priority[i]].qp = qps_;
      rnn_.insert(std::make_pair(qps_, arm_priority[i]));
      qps_++;
    }
  }

  for (int i=0; i < num_rel_type_; i++) {
    std::vector<RobotRelation> r_vec_;
    for (int j=0; j < num_rel_; j++) {
      RobotRelation rel_;
      rel_.base_arm = yn_["relation_types"][i]["relations"][j]["arm_set"][0].as<std::string>();
      rel_.target_arm = yn_["relation_types"][i]["relations"][j]["arm_set"][1].as<std::string>();
      rel_.relation = Eigen::Vector3d(yn_["relation_types"][i]["relations"][j]["relation"].as<std::vector<double>>().data());
      r_vec_.push_back(rel_);
    }
    file_names_.push_back(yn_["relation_types"][i]["file_name"].as<std::string>());
    rr_.push_back(r_vec_);
  }
  wfi_.open(dws_ + save_name_ + "_iteration_info.txt");
  midpoint_path = dws_ + save_name_ + "_midpoint_save.txt";
  wfi_ << "description: " <<  yn_["description"].as<std::string>() << std::endl;
  wfi_ << "calibrate base: " << calibrate_base_ << std::endl;
  wfi_ << "continue calib: " << continue_calib_ << std::endl;
  wfi_ << "method: " << method_ << std::endl;
  wfi_ << "lambda: " << lambda_ << std::endl;
  std::cout << "description: " <<  yn_["description"].as<std::string>() << std::endl;
  std::cout << "save name: " << save_name_ << std::endl;
  std::cout << "calibrate base: " << calibrate_base_ << std::endl;
  std::cout << "continue calib: " << continue_calib_ << std::endl;
  std::cout << "method: " << method_ << std::endl;
  std::cout << "lambda: " << lambda_ << std::endl;
}

void ClosedCalibration::structData() {
  panda_dh <<  // a       d     theta  alpha
                 0.0,    0.333,  0.0,  0.0,
                 0.0,    0.0,    0.0,  -M_PI_2,
                 0.0,    0.316,  0.0,  M_PI_2,
                 0.0825, 0.0,    0.0,  M_PI_2,
                 -0.0825,0.384,  0.0,  -M_PI_2,
                 0.0,    0.0,    0.0,  M_PI_2,
                 0.088,  0.0,    0.0,  M_PI_2;

  tot_j_ = num_arm_ * num_j_;
  input_size_stack_.push_back(0);
  for (int i=0; i< num_rel_type_; i++) {
    rf_.open(iws_ + file_names_[i]);
    std::vector<Eigen::VectorXd> inp_;
    while (!rf_.eof()) {
      Eigen::VectorXd d(tot_j_);
      for (int idx=0; idx < tot_j_; idx++) {
        rf_ >> d(idx);
      }
      inp_.push_back(d);
    }
    inputs.push_back(inp_);
    num_input_ += inp_.size();
    input_size_stack_.push_back(num_input_);
    rf_.close();
  }
  del_pos_.resize(num_rel_ * num_input_ * per_rel_);
  del_dh_.resize(num_arm_ * num_j_ * num_dh_);
  dh_mat_.resize(num_arm_ * num_j_, num_dh_);
  jacobian_.resize(num_rel_ * num_input_ * per_rel_, num_arm_ * num_j_ * num_dh_);
  del_pos_.setZero();
  del_dh_.setZero();
  dh_mat_.setZero();
  jacobian_.setZero();
}

void ClosedCalibration::loadData() {
    rf_.open(midpoint_path);
    for (int arm_=0; arm_ < num_arm_; arm_++) {
      for (int j=0; j < num_j_; j++) {
        for (int d=0; d < num_dh_; d++) {
          rf_ >> dh_mat_(j + arm_ * num_j_, d);
        }
      }
    }
    rf_.close();
}

Eigen::VectorXd ClosedCalibration::computeDeltaPos() {
  Eigen::VectorXd rtn_del_pos(num_rel_ * num_input_ * per_rel_);
  for (int ty=0; ty < num_rel_type_; ty++) {
    for (int iq=0; iq < inputs[ty].size(); iq++) {
      for (int ir=0; ir < num_rel_; ir++) {
        const std::string &b_ = rr_[ty][ir].base_arm;
        const std::string &t_ = rr_[ty][ir].target_arm;
        const Eigen::Ref<const Eigen::Vector3d> &tr_ = rr_[ty][ir].relation;  // true relation
        const Eigen::Ref<const Eigen::VectorXd> &bq_ = inputs[ty][iq].segment<num_j_>(num_j_ * rs_[b_].qp);
        const Eigen::Ref<const Eigen::VectorXd> &tq_ = inputs[ty][iq].segment<num_j_>(num_j_ * rs_[t_].qp);
        Eigen::Isometry3d bas_ =  rs_[b_].base * rs_[b_].fpm.getTransform(bq_);
        Eigen::Isometry3d tar_ =  rs_[t_].base * rs_[t_].fpm.getTransform(tq_);
        Eigen::Isometry3d T_bt_ = bas_.inverse() * tar_;
        int sp_ = ir*per_rel_ + iq*num_rel_*per_rel_ + input_size_stack_[ty]*num_rel_*per_rel_;
        rtn_del_pos.segment<per_rel_>(sp_) = tr_ - T_bt_.translation();
        #ifdef DEBUG_MODE
        // Eigen::Isometry3d test_bas_ = rs_[b_].base * computePandaTransform(panda_dh, dh_mat_.block<num_j_, num_dh_>(rs_[b_].qp * num_j_, 0), bq_);
        // Eigen::Isometry3d test_tar_ = rs_[t_].base * computePandaTransform(panda_dh, dh_mat_.block<num_j_, num_dh_>(rs_[t_].qp * num_j_, 0), tq_);
        // Eigen::Isometry3d test_T_ = test_bas_.inverse() * test_tar_;
        // std::cout << "\n\nDEBUG\n\n";
        // std::cout << "base_test:\n" << test_bas_.matrix() << std::endl;
        // std::cout << "base_modl:\n" << bas_.matrix() << std::endl;
        // std::cout << "target_test:\n" << test_tar_.matrix() << std::endl;
        // std::cout << "target_modl:\n" << tar_.matrix() << std::endl;
        // std::cout << "\ntranslation     : " << T_bt_.translation().transpose() << std::endl;
        // std::cout << "translation_test: " << test_T_.translation().transpose() << std::endl;
        // std::cout << "relation        : " << tr_.transpose() << std::endl;
        // std::cout << "\ntranslation.norm()     : " << T_bt_.translation().norm() << std::endl;
        // std::cout << "translation_test.norm(): " << test_T_.translation().norm() << std::endl;
        // std::cout << "relation.norm()        : " << tr_.norm() << std::endl;
        #endif
      }
    }
  }
  return rtn_del_pos;
}

void ClosedCalibration::updateStatus() {
  for (int ai=0; ai < num_arm_; ai++) {
    const Eigen::Ref<const Eigen::Matrix<double, num_j_, num_dh_>> &mat_ = dh_mat_.block<num_j_, num_dh_>(ai * num_j_, 0);
    std::cout << "\ninitializing " << rnn_[ai] << " with dh:\n" << mat_ << std::endl;
    rs_[rnn_[ai]].fpm.initModel(mat_);
  }
  std::cout << std::endl << std::endl;
}

void ClosedCalibration::computeJacobian() {
  // Use a 7-point central difference stencil method.
  Eigen::VectorXd t1, t2, m1, m2, m3;
  for (int arm_=0; arm_< num_arm_; arm_++) {
    Eigen::Matrix<double, num_j_, num_dh_> arm_dh_ = dh_mat_.block<num_j_, num_dh_>(arm_ * num_j_, 0);
    FrankaPandaModel &arm_fpm_ = rs_[rnn_[arm_]].fpm;
    Eigen::Matrix<double, num_j_, num_dh_> y1 = arm_dh_;
    Eigen::Matrix<double, num_j_, num_dh_> y2 = arm_dh_;
    for (int j=0; j< num_j_; j++) {
      for (int dh_=0; dh_< num_dh_; dh_++) {
        const double ax = std::fabs(arm_dh_(j, dh_));
        // Make step size as small as possible while still giving usable accuracy.
        const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);
        y1(j, dh_) += h;
        arm_fpm_.initModel(y1);
        t1 = computeDeltaPos();
        y2(j, dh_) -= h;
        arm_fpm_.initModel(y2);
        t2 = computeDeltaPos();
        m1 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        y1(j, dh_) += h;
        arm_fpm_.initModel(y1);
        t1 = computeDeltaPos();
        y2(j, dh_) -= h;
        arm_fpm_.initModel(y2);
        t2 = computeDeltaPos();
        m2 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        y1(j, dh_) += h;
        arm_fpm_.initModel(y1);
        t1 = computeDeltaPos();
        y2(j, dh_) -= h;
        arm_fpm_.initModel(y2);
        t2 = computeDeltaPos();
        m3 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        jacobian_.col(arm_*num_j_*num_dh_ + j*num_dh_ + dh_) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

        // Reset for next iteration.
        y1(j, dh_) = y2(j, dh_) = arm_dh_(j, dh_);
      }
    }

    // Reset for next iteration.
    arm_fpm_.initModel(arm_dh_);
  }
}

void ClosedCalibration::executeCalibration() {
  const int max_iter = 10000;
  int iter = max_iter;
  double ev_b;
  double ev_f = 0.0;
  double rate_current;
  int min_counter = 0;
  int inv_counter = 0;

  while (iter--) {
    std::cout << "\n\n----------------------------------\niteration: " << max_iter - iter << "\n" << std::endl;
    wfi_ << "\n\n----------------------------------\niteration: " << max_iter - iter << "\n" << std::endl;

    updateStatus();
    ev_b = ev_f;
    del_pos_ = computeDeltaPos();
    ev_f = del_pos_.squaredNorm() / num_input_;
    computeJacobian();
    if (iter < max_iter - 1) {
      rate_current = ((ev_b - ev_f) / ev_b) * 100.0;
      std::cout << "rate: " << rate_current << std::endl;
      wfi_ << "rate: " << rate_current << std::endl;
      if (0 < rate_current < 0.015) {
        min_counter++;
        if (min_counter == 2) {
          lambda_ /= lambda_factor_;
          std::cout << "Lambda(CHANGED): " << lambda_ << " ----------------------------------" << std::endl;
          wfi_      << "Lambda(CHANGED): " << lambda_ << " ----------------------------------" << std::endl;
          min_counter = 0;
        }
      } else if (rate_current < -0.01) {
        inv_counter++;
        if (inv_counter == 2) {
          lambda_ *= lambda_factor_;
          std::cout << "Lambda(CHANGED): " << lambda_ << " ----------------------------------" << std::endl;
          wfi_      << "Lambda(CHANGED): " << lambda_ << " ----------------------------------" << std::endl;
          inv_counter = 0;
        }
      }
    }

    std::cout << "eval(before): " << ev_f << std::endl;

    Eigen::MatrixXd weight(num_arm_ * num_j_ * num_dh_, num_arm_ * num_j_ * num_dh_);
    weight.setIdentity();
    const double rate_ = 1e-2;
    for (int ai=0; ai < num_arm_; ai++) {
      int arm_idx_ = ai * num_j_ * num_dh_;
      // weight(0 + arm_idx_, 0 + arm_idx_) = rate_;
      weight(1 + arm_idx_, 1 + arm_idx_) = rate_;
      weight(25 + arm_idx_, 25 + arm_idx_) = rate_;
    }

    if (method_ == "lm_method") {
      auto & j = jacobian_;
      Eigen::MatrixXd j_diag = (j.transpose() * j).diagonal().asDiagonal();
      Eigen::MatrixXd identity_ = j_diag;
      auto j_inv = (j.transpose() * j + lambda_ * (j_diag + identity_.setIdentity())).inverse() * j.transpose();
      del_dh_ = weight * j_inv * del_pos_;
      #ifdef DEBUG_MODE
      // std::cout << "jacobian_:\n" << jacobian_.matrix() << std::endl;
      // std::cout << "j_inv:\n" << j_inv.matrix() << std::endl;
      // std::cout << "j_diag  : " << j_diag.diagonal().transpose() << std::endl;
      // std::cout << "del_dh_ : " << del_dh_.transpose() << std::endl;
      #endif
    } else if (method_ == "svd_method") {
      // Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian_, Eigen::ComputeThinU | Eigen::ComputeThinV);
      // std::cout << "testing svd\nsingular values are:\n" << svd.singularValues() << std::endl;
      // del_dh_ = weight * svd.solve(del_pos_);
      del_dh_ = weight * jacobian_.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(del_pos_);
    }

    std::cout << "del_dh_.norm(): " << del_dh_.norm() << std::endl;

    for (int arm_=0; arm_< num_arm_; arm_++) {
      for (int j=0; j< num_j_; j++) {
        dh_mat_.row(arm_ * num_j_ + j).head<num_dh_>() -= del_dh_.segment<num_dh_>((arm_ * num_j_ + j)*num_dh_);
      }
    }

    std::ofstream offset_save(midpoint_path);
    wfi_ << "eval(before): "<< del_pos_.squaredNorm() / num_input_ << std::endl;
    wfi_ << "del_dh_.norm(): "<< del_dh_.norm() << std::endl;
    for (int ai=0; ai < num_arm_; ai++) {
      std::string nm_ = rs_[rnn_[ai]].name;
      Eigen::Matrix<double, num_j_, num_dh_> mat_ = dh_mat_.block<num_j_, num_dh_>(ai * num_j_, 0);
      wfi_ << nm_ << std::endl;
      wfi_ << mat_.format(tab_format) << std::endl;
      offset_save << mat_.format(tab_format);
      if (ai < num_arm_-1)
        offset_save << std::endl;
    }
    wfi_ <<"\n"<< std::endl;
    offset_save.close();
  }
}
