// Copyright 2021 Kim Hyoung Cheol. All rights reserved.
// Author: Kim Hyoung Cheol
// Email : kimhc37@snu.ac.kr

#include "../include/fixed_calibration.h"
// #define TEST_PRINT

const int N_DH = 4;  // number of dh parameters to calibrate
const int N_J = 7;  // number of joints to calibrate
const int N_CAL = N_DH*N_J;  // total variables to calibrate

std::vector<Eigen::Matrix<double, N_J, 1>> q_input_1, q_input_2, q_input_3;
Eigen::VectorXd dh_al(N_J), dh_a(N_J), dh_d(N_J);
Eigen::Vector3d true_p_1;
Eigen::Vector3d true_p_2;
Eigen::Vector3d true_p_3;
Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");

std::string current_workspace = "/home/kimhc/git/kinematics_calibration/data/fixed_calibration/";
std::string data_input;
std::string arm_name = "panda_left";
int fixed_points = 2;
std::ofstream iteration_info(current_workspace + "debug/" + arm_name + "_iteration_info.txt");
Eigen::Isometry3d dh_to_transform(const double a, const double d, const double alpha, const double theta)
{
  Eigen::Isometry3d transform_dh;
  transform_dh.setIdentity();
  transform_dh.linear()       << cos(theta),             -1*sin(theta),          0.0,
                                 sin(theta)*cos(alpha),  cos(theta)*cos(alpha),  -1*sin(alpha),
                                 sin(theta)*sin(alpha),  cos(theta)*sin(alpha),  cos(alpha);
  transform_dh.translation()  << a, -1*sin(alpha)*d, cos(alpha)*d;
  return transform_dh;
}

std::vector<Eigen::Matrix<double, N_J, 1>> read_data(std::string data_name)
{
  std::ifstream rf;
  std::vector<Eigen::Matrix<double, N_J, 1>> q_input_;
  data_input = current_workspace + "input_data/" + data_name + ".txt";
  std::cout << "reading -- " << data_input << std::endl; 
  rf.open(data_input);
  while (!rf.eof())
  {
    Eigen::Matrix<double, N_J, 1> d;
    for (int i=0; i<N_J; i++)
    {
      rf >> d(i);
    }
    q_input_.push_back(d);
  }
  rf.close();
  std::cout << "complete - size: " << q_input_.size() << std::endl;
  return q_input_;
}

void write_pos_info(std::string data_name, std::vector<Eigen::Matrix<double, N_J, 1>> q_input_, Eigen::Matrix<double, N_J, 4> dh)
{
  std::ofstream x_out(current_workspace + "debug/" + data_name + ".txt");
  auto fpm = FrankaPandaModel();
  fpm.initModel(dh);

  for (auto & q : q_input_)
  {
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    T_0e = fpm.getTransform(q);
    x_out << T_0e.translation().transpose() << std::endl;
  }
}

int main(int argc, char**argv)
{
  std::cout << "calibrating arm: " << arm_name << std::endl;
  iteration_info << "calibrating arm: " << arm_name << std::endl;
  // 0.07(from 7 to bottom) - 0.006(from robot base to bottom)
  double z = 0.064;
  if (arm_name == "panda_left")
  {
    q_input_1 = read_data("input_data_left_1");
    q_input_2 = read_data("input_data_left_2");
    if (fixed_points == 3) q_input_3 = read_data("input_data_left_3");
    true_p_1 << -0.075, -0.4, z;
    true_p_2 << -0.5, 0.0, z;
    if (fixed_points == 3) true_p_3 << 0.675, 0, z;
  }
  else if (arm_name == "panda_right")
  {
    q_input_1 = read_data("input_data_right_1");
    q_input_2 = read_data("input_data_right_2");
    if (fixed_points == 3) q_input_3 = read_data("input_data_right_3");
    true_p_1 << -0.075, 0.4, z;
    true_p_2 << -0.5, 0.0, z;
    if (fixed_points == 3) true_p_3 << 0.675, 0, z;
  }
  std::cout << "true position 1: " << true_p_1.transpose() << std::endl;
  iteration_info << "true position 1: " << true_p_1.transpose() << std::endl;
  std::cout << "true position 2: " << true_p_2.transpose() << std::endl;
  iteration_info << "true position 2: " << true_p_2.transpose() << std::endl;
  if (fixed_points == 3)
  {
  std::cout << "true position 3: " << true_p_3.transpose() << std::endl;
  iteration_info << "true position 3: " << true_p_3.transpose() << std::endl;
  }

  typedef std::pair < Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 3, 1> > caldata ;

  std::vector <caldata> calib_dataset;

  // memory prealloc
  if (fixed_points == 2) calib_dataset.reserve(q_input_1.size() + q_input_2.size());
  else if (fixed_points == 3) calib_dataset.reserve(q_input_1.size() + q_input_2.size() + q_input_3.size());

  for (auto & q : q_input_1)
  {
    calib_dataset.push_back(std::make_pair(q, true_p_1));
  }
  for (auto & q : q_input_2)
  {
    calib_dataset.push_back(std::make_pair(q, true_p_2));
  }
  if (fixed_points == 3)
  {
    for (auto & q : q_input_3)
    {
      calib_dataset.push_back(std::make_pair(q, true_p_3));
    }
  }

  dh_al << 0.0, -1.0*M_PI_2, M_PI_2, M_PI_2, -1.0*M_PI_2, M_PI_2, M_PI_2;
  dh_a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
  dh_d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;

  auto fpm = FrankaPandaModel();

  auto function_kim2 = [&calib_dataset,&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::VectorXd> out) {
    const auto & q = c.first;
    Eigen::Matrix<double, N_J, 4> dh;
    dh.setZero();
    for (int i=0 ; i<N_J; i++){dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);}

    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    fpm.initModel(dh);
    T_0e = fpm.getTransform(q);
    auto t = T_0e.translation();

    out = c.second - t;
  };

  auto jacobian_kim2 = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::MatrixXd> out)
  {
    const auto & q = c.first;
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie, p_ie_1;
    Eigen::Matrix3d R_0i, R_0i_1;
    Eigen::Isometry3d T_0i, T_0e, T_ie;
    Eigen::Matrix<double, 3, N_CAL> jacob_k; jacob_k.setZero();
    x_0i << 1.0, 0.0, 0.0; y_0i << 0.0, 1.0, 0.0; z_0i << 0.0, 0.0, 1.0;
    R_0i.setIdentity(); T_0e.setIdentity(); T_0i.setIdentity();
    Eigen::Matrix<double, N_J, 4> dh;
    dh.setZero();
    for (int i=0 ; i<N_J; i++){dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);}
    Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
    Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);

    fpm.initModel(dh);
    T_0e = fpm.getTransform(q);

    Eigen::Isometry3d T_0e_test;
    T_0e_test.setIdentity();
    for (int i=0; i<7; i++)
    {
      T_0e_test = T_0e_test * dh_to_transform(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i) + alpha_offset(i), q(i) + q_offset(i));
    }
    T_0e_test = T_0e_test * dh_to_transform(0.0, 0.107, 0.0, 0.0);

#ifdef TEST_PRINT
    std::cout << "--------------------------------------------------------------" << std::endl;
    std::cout << "true position: " << c.second.transpose() << std::endl;
    std::cout << "forward kinematics position: " << T_0e.translation().transpose() << std::endl;
    // std::cout << "franka model updater:\n" << T_0e.matrix() << std::endl;
    // std::cout << "function made:\n" << T_0e_test.matrix() << std::endl;
    std::cout << "translation diff: " << (T_0e.translation() - c.second).norm() << std::endl;

#endif

    p_ie = T_0e.translation();

    for (int i=0; i<N_J; i++)
    {
      x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i; p_ie_1 = p_ie; R_0i_1 = R_0i;

      x_0i = cos(q(i)+q_offset(i)) * x_0i_1               + sin(q(i)+q_offset(i)) * cos(dh_al(i)+alpha_offset(i)) * y_0i_1 + sin(q(i)+q_offset(i)) * sin(dh_al(i)+alpha_offset(i)) * z_0i_1;
      y_0i = -1.0*sin(q(i)+q_offset(i)) * x_0i_1          + cos(q(i)+q_offset(i)) * cos(dh_al(i)+alpha_offset(i)) * y_0i_1 + cos(q(i)+q_offset(i)) * sin(dh_al(i)+alpha_offset(i)) * z_0i_1;
      z_0i = -1.0*sin(dh_al(i)+alpha_offset(i)) * y_0i_1  + cos(dh_al(i)+alpha_offset(i)) * z_0i_1;

      T_0i = T_0i * dh_to_transform(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i), dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
      R_0i = T_0i.linear();
      p_ie = (T_0i.inverse() * T_0e).translation();

      jacob_k.col(0 + i*N_DH) += x_0i_1;
      jacob_k.col(1 + i*N_DH) += z_0i;
      jacob_k.col(2 + i*N_DH) += z_0i.cross(R_0i*p_ie);
      if (N_DH == 4){jacob_k.col(3 + i*N_DH) += x_0i_1.cross(R_0i_1*p_ie_1);}
    }
    out = -jacob_k;
  };
  Eigen::VectorXd x(N_CAL);
  x.setZero();

  double eval_current, eval_before, rate_current;
  int total_len = q_input_1.size() + q_input_2.size();
  if (fixed_points == 3) total_len += q_input_3.size();
  int max_iter = 10000;
  int iter = max_iter;
  int end_counter = 0;

  Eigen::VectorXd p_total;
  Eigen::MatrixXd jac_total;
  Eigen::VectorXd del_phi;
  p_total.resize(3*total_len);
  jac_total.resize(3*total_len, N_CAL);
  del_phi.resize(N_CAL);

  double lambda = 1.0;
  int min_counter = 0;

  while (iter--)
  {
    for (int i=0; i<calib_dataset.size(); ++i)
    {
      function_kim2(x, calib_dataset[i], p_total.segment<3>(i*3));
      jacobian_kim2(x, calib_dataset[i], jac_total.block<3,N_CAL>(i*3, 0));
    }

    eval_before = eval_current;
    eval_current = p_total.norm() / total_len;
    if (iter == max_iter - 1) eval_before = eval_current;
    rate_current = ((eval_before - eval_current) / eval_before) * 100.0;

    // del_phi = jac_total.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p_total);

    Eigen::Matrix<double, N_CAL, N_CAL> weight;
    weight.setIdentity();
    const double rate_ = 1e-2;
    // weight(25,25) = rate_;
    // weight(25 + N_JDH*1,25 + N_JDH*1) = rate_;
    // if (N_ARM == 3) weight(25 + N_JDH*2,25 + N_JDH*2) = rate_;
    // weight(1,1) = rate_;
    // weight(1 + N_JDH*1,1 + N_JDH*1) = rate_;

    auto & j = jac_total;

    // LM method
    Eigen::MatrixXd j_diag = (j.transpose() * j ).diagonal().asDiagonal();
    auto j_inv = (j.transpose() * j + lambda * j_diag).inverse() * j.transpose();
    del_phi = weight * j_inv * p_total;


    x -= del_phi; // jacobi is oppisite direction

    if (iter % 10 == 0)
    {
      std::cout << "\n----------------------------------------\niter: " << max_iter - iter << std::endl;
      std::cout << "eval: " << eval_current << std::endl;
      std::cout << "rate: " << rate_current << std::endl;
      std::cout << "dphi: " << del_phi.norm() << std::endl;
      iteration_info << "\n----------------------------------------\niter: " << max_iter - iter << std::endl;
      iteration_info << "eval: " << eval_current << std::endl;
      iteration_info << "rate: " << rate_current << std::endl;
      iteration_info << "dphi: " << del_phi.norm() << std::endl;
    }

    if (rate_current < 0.0003)
    {
      min_counter ++;
      if (min_counter == 2)
      {
        lambda += 1.0;
        std::cout << "Lambda(CHANGED): " << lambda << " ----------------------------------" << std::endl;
        iteration_info << "Lambda(CHANGED): " << lambda << " ----------------------------------" << std::endl;
        min_counter = 0;
      }
    }
  }

  Eigen::Matrix<double, N_J, 4> dh;
  dh.setZero();
  for (int i=0 ; i<N_J; i++)
  {
    dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);
  }

  write_pos_info("test1", q_input_1, dh);
  write_pos_info("test2", q_input_2, dh);
  if (fixed_points == 3) write_pos_info("test3", q_input_3, dh);

  std::ofstream dh_out(current_workspace + "result/" + arm_name + "_dh_output.txt");

  dh_out << dh.format(tab_format);

  return 0;
}