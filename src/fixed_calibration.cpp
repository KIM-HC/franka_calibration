// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#include "../include/fixed_calibration.h"

int main(int argc, char**argv) {
  dataSetStruct data_set;
  FrankaPandaModel fp_model = FrankaPandaModel();

  if (argc == 2) {
    data_set.yaml_path = std::string(ws_) + "yaml/" + argv[1] + ".yaml";
  } else {
    data_set.yaml_path = std::string(ws_) + "yaml/fixed_left_1.yaml";
  }

  readYaml(&data_set);
  structData(&data_set);

  double eval_current, eval_before, rate_current;
  int max_iter = 10000;
  int iter = max_iter;
  int end_counter = 0;
  int min_counter = 0;

  Eigen::Matrix<double, N_CAL, N_CAL> weight;
  weight.setIdentity();
  const double rate_ = 1e-2;
  weight(25,25) = rate_;
  weight(1,1) = rate_;

  while (iter--) {
    computeDeltaPos(&data_set, &fp_model);
    computeJacobian(&data_set, &fp_model);

    eval_before = eval_current;
    eval_current = data_set.del_pos.norm() / data_set.num_input;
    if (iter == max_iter - 1) eval_before = eval_current;
    rate_current = ((eval_before - eval_current) / eval_before) * 100.0;

    if (data_set.method == "lm_method") {
      auto & j = data_set.jacobian;
      Eigen::MatrixXd j_diag = (j.transpose() * j).diagonal().asDiagonal();
      auto j_inv = (j.transpose() * j + lambda * j_diag).inverse() * j.transpose();
      data_set.del_dh = weight * j_inv * data_set.del_pos;
      if (rate_current < 0.0003) {
        min_counter++;
        if (min_counter == 2) {
          lambda *= 1.05;
          std::cout << "Lambda(CHANGED): " << lambda << " ----------------------------------" << std::endl;
          iteration_info << "Lambda(CHANGED): " << lambda << " ----------------------------------" << std::endl;
          min_counter = 0;
        }
      }
    } else if (data_set.method == "svd_method") {
      data_set.del_dh = weight * data_set.jacobian.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(data_set.del_pos);
    }

    data_set.dh_vec -= data_set.del_dh;  // jacobi is oppisite direction
    for (int i=0 ; i< N_J; i++)
      data_set.dh_mat.row(i).head<N_DH>() = data_set.dh_vec.segment<N_DH>(i*N_DH);
    fp_model.initModel(data_set.dh_mat);

    if ((max_iter - iter) < 20) {
      mid_point_save << "\n----------------------------------------\niter: " << max_iter - iter << std::endl;
      mid_point_save << "eval: " << eval_current << std::endl;
      mid_point_save << "rate: " << rate_current << std::endl;
      mid_point_save << "del_dh: " << data_set.del_dh.norm() << std::endl;
      mid_point_save << data_set.dh_mat.format(tab_format) << std::endl;
    }
    if (iter % 10 == 0) {
      std::cout << "\n----------------------------------------\niter: " << max_iter - iter << std::endl;
      std::cout << "eval: " << eval_current << std::endl;
      std::cout << "rate: " << rate_current << std::endl;
      std::cout << "del_dh: " << data_set.del_dh.norm() << std::endl;
      iteration_info << "\n----------------------------------------\niter: " << max_iter - iter << std::endl;
      iteration_info << "eval: " << eval_current << std::endl;
      iteration_info << "rate: " << rate_current << std::endl;
      iteration_info << "del_dh: " << data_set.del_dh.norm() << std::endl;
      if ((iter % 30 == 0) && ((max_iter - iter) > 20)) {
        mid_point_save << "\n----------------------------------------\niter: " << max_iter - iter << std::endl;
        mid_point_save << "eval: " << eval_current << std::endl;
        mid_point_save << "rate: " << rate_current << std::endl;
        mid_point_save << "del_dh: " << data_set.del_dh.norm() << std::endl;
        mid_point_save << data_set.dh_mat.format(tab_format) << std::endl;
      }
    }
  }

  writePosInfo(&data_set, &fp_model);

  std::ofstream ofs(std::string(ws_) + "result/" + data_set.arm_name
                    + "_" + data_set.method + "_" + std::to_string(sn) + "_dh_output.txt");
  ofs << data_set.dh_mat.format(tab_format);

  ofs.close();
  iteration_info.close();
  mid_point_save.close();
  return 0;
}

void computeJacobian(dataSetStruct *ds, FrankaPandaModel *fpm) {
  Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie, p_ie_1;
  Eigen::Matrix3d R_0i, R_0i_1;
  Eigen::Isometry3d T_0i, T_0e, T_ie;
  Eigen::Matrix<double, 3, N_CAL> jacob_k;
  Eigen::Ref<Eigen::VectorXd> a_offset = ds->dh_mat.col(0);
  Eigen::Ref<Eigen::VectorXd> d_offset = ds->dh_mat.col(1);
  Eigen::Ref<Eigen::VectorXd> q_offset = ds->dh_mat.col(2);
  Eigen::Ref<Eigen::VectorXd> alpha_offset = ds->dh_mat.col(3);

  for (int i=0; i< ds->inputs.size(); i++) {
    jacob_k.setZero();
    x_0i << 1.0, 0.0, 0.0;
    y_0i << 0.0, 1.0, 0.0;
    z_0i << 0.0, 0.0, 1.0;
    R_0i.setIdentity();
    T_0e.setIdentity();
    T_0i.setIdentity();

    Eigen::Matrix<double, N_J, 1> &q = ds->inputs[i].first;
    T_0e = fpm->getTransform(q);
    p_ie = T_0e.translation();

    for (int j=0; j< N_J; j++) {
      x_0i_1 = x_0i;
      y_0i_1 = y_0i;
      z_0i_1 = z_0i;
      p_ie_1 = p_ie;
      R_0i_1 = R_0i;

      x_0i = cos(q(j)+q_offset(j)) * x_0i_1
              + sin(q(j)+q_offset(j)) * cos(panda_dh(j, 3)
              + alpha_offset(j)) * y_0i_1
              + sin(q(j)+q_offset(j)) * sin(panda_dh(j, 3)+alpha_offset(j)) * z_0i_1;

      y_0i = -1.0*sin(q(j)+q_offset(j)) * x_0i_1
              + cos(q(j)+q_offset(j)) * cos(panda_dh(j, 3)
              + alpha_offset(j)) * y_0i_1
              + cos(q(j)+q_offset(j)) * sin(panda_dh(j, 3)+alpha_offset(j)) * z_0i_1;

      z_0i = -1.0*sin(panda_dh(j, 3)+alpha_offset(j)) * y_0i_1
              + cos(panda_dh(j, 3)+alpha_offset(j)) * z_0i_1;

      T_0i = T_0i * transformDH(panda_dh(j, 0) + a_offset(j),
                               panda_dh(j, 1) + d_offset(j),
                               panda_dh(j, 3) + alpha_offset(j),
                               q(j) + q_offset(j));
      R_0i = T_0i.linear();
      p_ie = (T_0i.inverse() * T_0e).translation();

      jacob_k.col(0 + j*N_DH) += x_0i_1;
      jacob_k.col(1 + j*N_DH) += z_0i;
      jacob_k.col(2 + j*N_DH) += z_0i.cross(R_0i*p_ie);
      if (N_DH == 4) {jacob_k.col(3 + j*N_DH) += x_0i_1.cross(R_0i_1*p_ie_1);}
    }
    ds->jacobian.block<3, N_CAL>(i*3, 0) = -jacob_k;
  }
}

void computeDeltaPos(dataSetStruct *ds, FrankaPandaModel *fpm) {
  for (int i=0; i< ds->inputs.size(); i++)
    ds->del_pos.segment<3>(i*3) = ds->inputs[i].second - fpm->getTransform(ds->inputs[i].first).translation();
}

void readYaml(dataSetStruct *ds) {
  YAML::Node yaml_reader = YAML::LoadFile(ds->yaml_path);
  ds->arm_name = yaml_reader["name"].as<std::string>();
  ds->method = yaml_reader["method"].as<std::string>();
  lambda = yaml_reader["lambda"].as<double>();
  sn = yaml_reader["save_number"].as<int>();
  ds->calibrate_base = yaml_reader["calibrate_base"];
  ds->num_ref = yaml_reader["relations"].size();

  iteration_info.open(std::string(ws_) + "debug/" + ds->arm_name
                      + "_" + std::to_string(sn) + "_iteration_info.txt");

  mid_point_save.open(std::string(ws_) + "debug/" + ds->arm_name
                      + "_" + std::to_string(sn) + "_dh_info.txt");

  iteration_info << "calibrate base: " << std::boolalpha << ds->calibrate_base << std::endl;
  iteration_info << "calibrating arm: " << ds->arm_name << std::endl;
  iteration_info << "method: " << ds->method << std::endl;
  std::cout << "calibrate base: " << std::boolalpha << ds->calibrate_base << std::endl;
  std::cout << "calibrating arm: " << ds->arm_name << std::endl;
  std::cout << "method: " << ds->method << std::endl;
  for (int i=0; i< ds->num_ref; i++) {
    ds->relations.push_back(Eigen::Vector3d(yaml_reader["relations"][i]["relation"].as<std::vector<double>>().data()));
    ds->file_names.push_back(yaml_reader["relations"][i]["file_name"].as<std::string>());
    iteration_info << "relation: " << i << ": " << ds->relations[i].transpose() << std::endl;
    iteration_info << "rel-path: " << i << ": " << ds->file_names[i] << std::endl;
    std::cout << "relation: " << i << ": " << ds->relations[i].transpose() << std::endl;
    std::cout << "rel-path: " << i << ": " << ds->file_names[i] << std::endl;
  }
}

void structData(dataSetStruct *ds) {
  panda_dh <<  // a       d     theta  alpha
                 0.0,    0.333,  0.0,  0.0,
                 0.0,    0.0,    0.0,  -M_PI_2,
                 0.0,    0.316,  0.0,  M_PI_2,
                 0.0825, 0.0,    0.0,  M_PI_2,
                 -0.0825,0.384,  0.0,  -M_PI_2,
                 0.0,    0.0,    0.0,  M_PI_2,
                 0.088,  0.0,    0.0,  M_PI_2;

  ds->base.setIdentity();
  ds->dh_mat.setZero();
  ds->dh_vec.setZero();

  for (int i=0; i< ds->num_ref; i++) {
    std::string data_input = std::string(ws_) + "input_data/" + ds->file_names[i];
    std::ifstream rf(data_input);
    while (!rf.eof()) {
      Eigen::Matrix<double, N_J, 1> d;
      for (int idx=0; idx < N_J; idx++) {
        rf >> d(idx);
      }
      ds->inputs.push_back(std::make_pair(d, ds->relations[i]));
    }
    rf.close();
  }
  ds->num_input = ds->inputs.size();
  ds->jacobian.resize(3*ds->num_input, N_CAL);
  ds->del_pos.resize(3*ds->num_input);
  std::cout << "num data: " << ds->num_input << std::endl;
}

void writePosInfo(dataSetStruct *ds, FrankaPandaModel *fpm) {
  std::vector<std::ofstream> writer;
  for (int i=0; i< ds->num_ref; i++) {
    writer.push_back(std::ofstream(std::string(ws_) + "debug/pos_info_" + "_"
                                   + ds->method + "_" + std::to_string(sn) + std::to_string(i) + ".txt"));
  }
  for (const auto &pair_ : ds->inputs) {
    Eigen::Isometry3d T_0e = fpm->getTransform(pair_.first);
    for (int i=0; i< ds->num_ref; i++) {
      if (pair_.second == ds->relations[i]) {
        writer[i] << T_0e.translation().transpose() << std::endl;
        break;
      }
    }
  }
}
