// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#include "../../include/utils/calibration_math.h"

namespace calib_math {
  Eigen::Isometry3d transformDH(const double &a,
                                    const double &d,
                                    const double &alpha,
                                    const double &theta) {
    Eigen::Isometry3d transform_dh_;
    transform_dh_.setIdentity();
    transform_dh_.linear()       << cos(theta),             -1*sin(theta),          0.0,
                                  sin(theta)*cos(alpha),  cos(theta)*cos(alpha),  -1*sin(alpha),
                                  sin(theta)*sin(alpha),  cos(theta)*sin(alpha),  cos(alpha);
    transform_dh_.translation()  << a, -1*sin(alpha)*d, cos(alpha)*d;
    return transform_dh_;
  }

  Eigen::Isometry3d computePandaTransform(const Eigen::Ref<const Eigen::Matrix<double, 7, 4>> dh,
                                    const Eigen::Ref<const Eigen::Matrix<double, 7, 4>> del_dh,
                                    const Eigen::Ref<const Eigen::Matrix<double, 7, 1>> q) {
    Eigen::Isometry3d pd_;
    pd_.setIdentity();
    for (int i=0; i<7; i++)
    {
      pd_ = pd_ * transformDH(dh(i, 0) + del_dh(i, 0),
                             dh(i, 1) + del_dh(i, 1),
                             dh(i, 3) + del_dh(i, 3),
                             q(i)     + del_dh(i, 2));
    }
    pd_ = pd_ * transformDH(0.0, 0.107, 0.0, -M_PI_4);
    return pd_;
  }
}  // namespace calib_math
