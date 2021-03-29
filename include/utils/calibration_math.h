// Copyright 2021 Kim Hyoung Cheol (kimhc37@snu.ac.kr). All rights reserved.

#ifndef UTILS_CALIBRATION_MATH_H_
#define UTILS_CALIBRATION_MATH_H_

#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <map>


namespace calib_math {

// receives dh parameters & returns homogeneous transform matrix
Eigen::Isometry3d dh_to_transform(const double a,
                                  const double d,
                                  const double alpha,
                                  const double theta);



}  // namespace calib_math

#endif  // UTILS_CALIBRATION_MATH_H_
