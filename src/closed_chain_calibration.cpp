// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#include "../include/closed_chain_calibration.h"

int main(int argc, char**argv) {
  if (argc >= 2) {
    std::string prefix;
    prefix = argv[1];
    data_iter = ws_ + data_ws_ + "debug/iteration_info_trial_" + prefix + ".txt";
    data_offset = ws_ + data_ws_ + "result/offset_data_trial_" + prefix + ".txt";
  } else {
    data_iter = ws_ + data_ws_ + "debug/iteration_info_trial_1.txt";
    data_offset = ws_ + data_ws_ + "result/offset_data_trial_1.txt";
  }
  std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Input path: " << ws_ + data_ws_ + "input_data/input_data.txt" << std::endl;
  std::cout << "Iter path: " << data_iter << std::endl;
  std::cout << "Offset path: " << data_offset << std::endl;
  std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
  initialize();

  if (argc >= 3) {
    rf.open(data_offset);
    for (int arm_=0; arm_< N_ARM; arm_++) {
      for (int j=0; j< N_J; j++) {
        for (int d=0; d< N_DH; d++) {
          rf >> offset_matrix[arm_](j, d);
        }
      }
    }
    rf.close();
    lambda = strtod(argv[2], NULL);
    std::cout << "Lambda: " << lambda << std::endl;
    std::cout << "\n<<USING SAVED OFFSET>>" << std::endl;
    std::cout << "offset_matrix LEFT:\n" << offset_matrix[0] << std::endl;
    std::cout << "\noffset_matrix RIGHT:\n" << offset_matrix[1] << std::endl;
    std::cout << "\noffset_matrix TOP:\n" << offset_matrix[2] << "\n\n" <<std::endl;
  } else {
    std::cout << "Lambda: " << lambda << std::endl;
  }

  std::ofstream iter_save(data_iter);
  const int max_iter = 10000;
  int iter = max_iter;
  double ev_b;
  double ev_f = 0.0;
  double rate_current;
  int min_counter = 0;

  while (iter--) {
    iter_save << "iteration: "<< max_iter - iter << std::endl;

    std::cout << "\n\n----------------------------------\niteration: " << max_iter - iter << "\n" << std::endl;
    for (int i=0; i< N_ARM; i++) {fpm[i].initModel(offset_matrix[i]);}
    ev_b = ev_f;
    p_total = getDistanceDiff();
    ev_f = p_total.squaredNorm() / num_data;
    getJacobian();
    if (iter < max_iter - 1) {
      rate_current = ((ev_b - ev_f) / ev_b) * 100.0;
      std::cout << "rate: " << rate_current << std::endl;
      iter_save << "rate: " << rate_current << std::endl;
      if (rate_current < 0.015) {
        min_counter++;
        if (min_counter == 2) {
          lambda += 1.0;
          std::cout << "Lambda(CHANGED): " << lambda << " ----------------------------------" << std::endl;
          iter_save << "Lambda(CHANGED): " << lambda << " ----------------------------------" << std::endl;
          min_counter = 0;
        }
      }
    }

    std::cout << "eval(before): " << ev_f << std::endl;

    Eigen::Matrix<double, N_CAL, N_CAL> weight;
    weight.setIdentity();
    const double rate_ = 1e-2;

    // weight(25,25) = rate_;
    // weight(25 + N_JDH*1,25 + N_JDH*1) = rate_;
    // if (N_ARM == 3) weight(25 + N_JDH*2,25 + N_JDH*2) = rate_;

    weight(1, 1) = rate_;
    weight(1 + N_JDH*1, 1 + N_JDH*1) = rate_;
    if (N_ARM == 3) weight(1 + N_JDH*2, 1 + N_JDH*2) = rate_;

    auto & j = jacobian;

    //// LM method
    Eigen::MatrixXd j_diag = (j.transpose() * j).diagonal().asDiagonal();
    auto j_inv = (j.transpose() * j + lambda * j_diag).inverse() * j.transpose();
    del_phi = weight * j_inv * p_total;
    // std::cout << "j_inv.diag:\n" << j_inv.diagonal() << std::endl;
    // std::cout << "j_diag.diag:\n" << j_diag.diagonal() << std::endl;

    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // std::cout << "testing svd\nsingular values are:\n" << svd.singularValues() << std::endl;
    // del_phi = svd.solve(p_total);
    // del_phi = weight * del_phi;

    std::cout << "del_phi.norm(): " << del_phi.norm() << std::endl;

    for (int arm_=0; arm_< N_ARM; arm_++) {
      for (int j=0; j< N_J; j++) {
        offset_matrix[arm_].row(j) -= del_phi.segment<N_DH>(arm_*N_JDH + j*N_DH);
      }
    }
#ifdef DEBUG_MODE
    std::cout << "\noffset_matrix LEFT:\n" << offset_matrix[0] << std::endl;
    std::cout << "\noffset_matrix RIGHT:\n" << offset_matrix[1] << std::endl;
    std::cout << "\noffset_matrix TOP:\n" << offset_matrix[2] << std::endl;
#endif
    // if (del_phi.norm() < 1e-9)
    // {
    //   std::cout<<"reached optimal value at iter: "<<100 - iter<<std::endl;
    //   break;
    // }

    iter_save << "eval(before): "<< p_total.squaredNorm() / num_data << std::endl;
    iter_save << "del_phi.norm(): "<< del_phi.norm() << std::endl;
    iter_save << "PANDA LEFT"<< std::endl;
    iter_save << offset_matrix[0].format(tab_format) << std::endl;
    iter_save << "PANDA RIGHT"<< std::endl;
    iter_save << offset_matrix[1].format(tab_format) << std::endl;
    iter_save << "PANDA TOP"<< std::endl;
    iter_save << offset_matrix[2].format(tab_format) << std::endl;
    iter_save <<"\n"<< std::endl;

    std::ofstream offset_save(data_offset);
    offset_save << offset_matrix[0].format(tab_format) <<std::endl;
    offset_save << offset_matrix[1].format(tab_format) <<std::endl;
    offset_save << offset_matrix[2].format(tab_format);
    offset_save.close();
  }

  return 0;
}

void read_data(std::string file_name) {
  rf.open(file_name);
  while (!rf.eof()) {
    Eigen::Matrix<double, N_ARM, N_J> ar;
    for (int i=0; i< N_ARM; i++) {
      Eigen::Matrix<double, 1, N_J> d;
      for (int j=0; j< N_J; j++) {
        rf >> d(j);
      }
      ar.row(i) = d;
    }
    int dn = 0;
    rf >> dn;
    data_number.push_back(dn);
    theta_data.push_back(ar);
  }
  rf.close();
}

void initialize() {
  read_data(ws_ + data_ws_ + "input_data/input_data.txt");

  num_data = theta_data.size();
  jacobian.resize(N_PARAM*num_data, N_CAL);
  p_total.resize(N_PARAM*num_data);

  z_rot180 << -1,  0,   0,
              0,   -1,  0,
              0,   0,   1;

  x_rot180 << 1,   0,   0,
              0,   -1,  0,
              0,   0,   -1;

  T_W0_LEFT.setIdentity();
  T_W0_LEFT.translation() << 0.0, 0.3, 1.0;
  T_W0_RIGHT.setIdentity();
  T_W0_RIGHT.translation() << 0.0, -0.3, 1.0;
  T_W0.push_back(T_W0_LEFT);
  T_W0.push_back(T_W0_RIGHT);
  if (N_ARM == 3) {
  T_W0_TOP.setIdentity();
  T_W0_TOP.translation() << 1.35, 0.3, 1.0;
  T_W0_TOP.linear() = T_W0_TOP.linear() * z_rot180;
  T_W0.push_back(T_W0_TOP);
  }

  std::vector<double> dist1;
  dist1.push_back(0.207);  // LEFT&RIGHT
  if (N_ARM == 3) {
    dist1.push_back(sqrt(pow(0.207, 2)*2));  // RIGHT&TOP
    dist1.push_back(0.207);  // TOP&LEFT
  }
  std::vector<double> dist2;
  dist2.push_back(sqrt(pow(0.207, 2) + pow(0.0051, 2)));  // LEFT&RIGHT
  if (N_ARM == 3) {
    dist2.push_back(sqrt(pow(0.207, 2)*2 + pow(0.0051, 2)));  // RIGHT&TOP
    dist2.push_back(0.207);  // TOP&LEFT
  }

  distTrue.push_back(dist1);
  distTrue.push_back(dist2);

  for (int i=0; i< N_ARM; i++) {
    FrankaPandaModel fpm_ = FrankaPandaModel();
    Eigen::Matrix<double, N_J, N_DH> x;
    x.setZero();

    offset_matrix.push_back(x);
    fpm.push_back(fpm_);
  }

  std::cout << "Number of datasets: " << theta_data.size() << std::endl;
  if (N_ARM == 2) {
    std::cout << "Distance L&R (data1): " << distTrue[0][0] << std::endl;
    std::cout << "Distance L&R (data2): " << distTrue[1][0] << std::endl;
  } else if (N_ARM == 3) {
    std::cout << "Distance L&R (data1): " << distTrue[0][0] << std::endl;
    std::cout << "Distance R&T (data1): " << distTrue[0][1] << std::endl;
    std::cout << "Distance T&L (data1): " << distTrue[0][2] << std::endl;
    std::cout << "Distance L&R (data2): " << distTrue[1][0] << std::endl;
    std::cout << "Distance R&T (data2): " << distTrue[1][1] << std::endl;
    std::cout << "Distance T&L (data2): " << distTrue[1][2] << std::endl;
  }
}

Eigen::VectorXd getDistanceDiff() {
  Eigen::VectorXd del_(N_PARAM*num_data);
  del_.setZero();

  for (int i=0; i< num_data; i++) {
    auto T_left = (T_W0[0]) * fpm[0].getTransform(theta_data[i].row(0));  // LEFT
    auto T_right = (T_W0[1]) * fpm[1].getTransform(theta_data[i].row(1));  // RIGHT
    Eigen::Isometry3d T_top;
    if (N_ARM == 3) T_top = (T_W0[2]) * fpm[2].getTransform(theta_data[i].row(2));  // TOP

    Eigen::Vector3d X_left = T_left.translation();
    Eigen::Vector3d X_right = T_right.translation();
    Eigen::Vector3d X_top;
    if (N_ARM == 3) X_top = T_top.translation();

    int nn = data_number[i] - 1;
    del_(i * N_PARAM + 0) = distTrue[nn][0] - (X_left - X_right).norm();
    if (N_ARM == 3) {
      del_(i * N_PARAM + 1) = distTrue[nn][1] - (X_right - X_top).norm();
      del_(i * N_PARAM + 2) = distTrue[nn][2] - (X_top - X_left).norm();
    }

    if (N_PARAM > N_ARM * (N_ARM-1) * 0.51) {
      if (nn == 1) {
        Eigen::Quaterniond q_left(T_left.linear());
        Eigen::Quaterniond q_right(T_right.linear());
        del_((i+0.5) * N_PARAM + 0) = M_PI - q_left.angularDistance(q_right);
        if (N_ARM == 3) {
          Eigen::Quaterniond q_top(T_top.linear() * z_rot180);
          del_((i+0.5) * N_PARAM + 1) = M_PI - q_right.angularDistance(q_top);
          del_((i+0.5) * N_PARAM + 2) = q_top.angularDistance(q_left);
        }
      } else {
        Eigen::Quaterniond q_left(T_left.linear());
        Eigen::Quaterniond q_right(T_right.linear());
        del_((i+0.5) * N_PARAM + 0) = q_left.angularDistance(q_right);
        if (N_ARM == 3) {
          Eigen::Quaterniond q_top(T_top.linear() * z_rot180);
          del_((i+0.5) * N_PARAM + 1) = q_right.angularDistance(q_top);
          del_((i+0.5) * N_PARAM + 2) = q_top.angularDistance(q_left);
        }
      }
    }

#ifdef DEBUG_MODE
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Data Number: " << data_number[i] << std::endl;
    std::cout << "Transform_LEFT: " << X_left.transpose() << std::endl;
    std::cout << "Transform_RIGHT: " << X_right.transpose() << std::endl;
    if (N_ARM == 3) {
      std::cout << "Transform_TOP: " << X_top.transpose() << std::endl;
    }
    std::cout << "(X_left - X_right).norm(): " << (X_left - X_right).norm() << std::endl;
    if (N_ARM == 3) {
      std::cout << "(X_right - X_top).norm(): " << (X_right - X_top).norm() << std::endl;
      std::cout << "(X_top - X_left).norm(): " << (X_top - X_left).norm() << std::endl;
    }
    std::cout << "del_(i * N_PARAM + 0): " << del_(i*N_PARAM + 0) << std::endl;
    if (N_ARM == 3) {
      std::cout << "del_(i * N_PARAM + 1): " << del_(i*N_PARAM + 1) << std::endl;
      std::cout << "del_(i * N_PARAM + 2): " << del_(i*N_PARAM + 2) << std::endl;
    }
    if (N_PARAM > N_ARM * (N_ARM-1) * 0.51) {
      std::cout << "del_((i+0.5) * N_PARAM + 0): " << del_((i+0.5) * N_PARAM + 0) << std::endl;
      if (N_ARM == 3) {
        std::cout << "del_((i+0.5) * N_PARAM + 1): " << del_((i+0.5) * N_PARAM + 1) << std::endl;
        std::cout << "del_((i+0.5) * N_PARAM + 2): " << del_((i+0.5) * N_PARAM + 2) << std::endl;
      }
    }
#endif
  }
  return del_;
}

void getJacobian() {
  // Use a 7-point central difference stencil method.
  Eigen::VectorXd t1, t2, m1, m2, m3;
  for (int arm_=0; arm_< N_ARM; arm_++) {
    Eigen::Matrix<double, N_J, N_DH> y1 = offset_matrix[arm_];
    Eigen::Matrix<double, N_J, N_DH> y2 = offset_matrix[arm_];
    for (int j=0; j< N_J; j++) {
      for (int dh_=0; dh_< N_DH; dh_++) {
        const double ax = std::fabs(offset_matrix[arm_](j, dh_));
        // Make step size as small as possible while still giving usable accuracy.
        const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);
        y1(j, dh_) += h;
        fpm[arm_].initModel(y1);
        t1 = getDistanceDiff();
        y2(j, dh_) -= h;
        fpm[arm_].initModel(y2);
        t2 = getDistanceDiff();
        m1 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        y1(j, dh_) += h;
        fpm[arm_].initModel(y1);
        t1 = getDistanceDiff();
        y2(j, dh_) -= h;
        fpm[arm_].initModel(y2);
        t2 = getDistanceDiff();
        m2 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        y1(j, dh_) += h;
        fpm[arm_].initModel(y1);
        t1 = getDistanceDiff();
        y2(j, dh_) -= h;
        fpm[arm_].initModel(y2);
        t2 = getDistanceDiff();
        m3 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        jacobian.col(arm_*N_JDH + j*N_DH + dh_) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

        // Reset for next iteration.
        y1(j, dh_) = y2(j, dh_) = offset_matrix[arm_](j, dh_);
      }
    }

    // Reset for next iteration.
    fpm[arm_].initModel(offset_matrix[arm_]);
  }
}
