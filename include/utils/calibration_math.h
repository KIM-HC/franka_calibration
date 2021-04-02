// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#pragma once
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <map>
#include <memory>  // std::shared_ptr

namespace calib_math {
    // receives dh parameters & returns homogeneous transform matrix
    Eigen::Isometry3d transformDH(const double &a,
                                 const double &d,
                                 const double &alpha,
                                 const double &theta);

    // computes transform of Franka Panda
    Eigen::Isometry3d computePandaTransform(const Eigen::Ref<const Eigen::Matrix<double, 7, 4>> dh,
                                            const Eigen::Ref<const Eigen::Matrix<double, 7, 4>> del_dh,
                                            const Eigen::Ref<const Eigen::Matrix<double, 7, 1>> q);

}  // namespace calib_math
