/*
 * interpolated_decay_grid.h
 *
 *  Created on: Sep 21, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_INTERPOLATED_DECAY_GRID_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_INTERPOLATED_DECAY_GRID_H_
#include <cmath>

#include "cartographer/mapping_3d/hybrid_grid.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {

// Interpolates between HybridGrid probability voxels. We use the tricubic
// interpolation which interpolates the values and has vanishing derivative at
// these points.
//
// This class is templated to work with the autodiff that Ceres provides.
// For this reason, it is also important that the interpolation scheme be
// continuously differentiable.
class InterpolatedDecayGrid {
 public:
  explicit InterpolatedDecayGrid(const HybridDecayGrid& hybrid_grid)
      : hybrid_grid_(hybrid_grid) {}

  InterpolatedDecayGrid(const InterpolatedDecayGrid&) = delete;
  InterpolatedDecayGrid& operator=(const InterpolatedDecayGrid&) = delete;

  // Returns the interpolated probability at (x, y, z) of the HybridGrid
  // used to perform the interpolation.
  //
  // This is a piecewise, continuously differentiable function. We use the
  // scalar part of Jet parameters to select our interval below. It is the
  // tensor product volume of piecewise cubic polynomials that interpolate
  // the values, and have vanishing derivative at the interval boundaries.
  double GetDecayRate(const double& x, const double& y, const double& z) const {
    double x1, y1, z1, x2, y2, z2;
    ComputeInterpolationDataPoints(x, y, z, &x1, &y1, &z1, &x2, &y2, &z2);

    const Eigen::Array3i index1(static_cast<int>(x), static_cast<int>(y),
                                static_cast<int>(z));
    const double q111 = hybrid_grid_.GetDecayRate(index1);
    const double q112 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(0, 0, 1));
    const double q121 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(0, 1, 0));
    const double q122 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(0, 1, 1));
    const double q211 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(1, 0, 0));
    const double q212 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(1, 0, 1));
    const double q221 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(1, 1, 0));
    const double q222 =
        hybrid_grid_.GetDecayRate(index1 + Eigen::Array3i(1, 1, 1));

    const double normalized_x = (x - x1) / (x2 - x1);
    const double normalized_y = (y - y1) / (y2 - y1);
    const double normalized_z = (z - z1) / (z2 - z1);

    // Compute pow(..., 2) and pow(..., 3). Using pow() here is very expensive.
    const double normalized_xx = normalized_x * normalized_x;
    const double normalized_xxx = normalized_x * normalized_xx;
    const double normalized_yy = normalized_y * normalized_y;
    const double normalized_yyy = normalized_y * normalized_yy;
    const double normalized_zz = normalized_z * normalized_z;
    const double normalized_zzz = normalized_z * normalized_zz;

    // We first interpolate in z, then y, then x. All 7 times this uses the same
    // scheme: A * (2t^3 - 3t^2 + 1) + B * (-2t^3 + 3t^2).
    // The first polynomial is 1 at t=0, 0 at t=1, the second polynomial is 0
    // at t=0, 1 at t=1. Both polynomials have derivative zero at t=0 and t=1.
    const double q11 = (q111 - q112) * normalized_zzz * 2. +
                  (q112 - q111) * normalized_zz * 3. + q111;
    const double q12 = (q121 - q122) * normalized_zzz * 2. +
                  (q122 - q121) * normalized_zz * 3. + q121;
    const double q21 = (q211 - q212) * normalized_zzz * 2. +
                  (q212 - q211) * normalized_zz * 3. + q211;
    const double q22 = (q221 - q222) * normalized_zzz * 2. +
                  (q222 - q221) * normalized_zz * 3. + q221;
    const double q1 = (q11 - q12) * normalized_yyy * 2. +
                 (q12 - q11) * normalized_yy * 3. + q11;
    const double q2 = (q21 - q22) * normalized_yyy * 2. +
                 (q22 - q21) * normalized_yy * 3. + q21;
    return (q1 - q2) * normalized_xxx * 2. + (q2 - q1) * normalized_xx * 3. +
           q1;
  }

 private:
  template <typename T>
  void ComputeInterpolationDataPoints(const T& x, const T& y, const T& z,
                                      double* x1, double* y1, double* z1,
                                      double* x2, double* y2,
                                      double* z2) const {
    const Eigen::Vector3f lower = CenterOfLowerVoxel(x, y, z);
    *x1 = lower.x();
    *y1 = lower.y();
    *z1 = lower.z();
    *x2 = lower.x() + 1;
    *y2 = lower.y() + 1;
    *z2 = lower.z() + 1;
  }

  // Center of the next lower voxel, i.e., not necessarily the voxel containing
  // (x, y, z). For each dimension, the largest voxel index so that the
  // corresponding center is at most the given coordinate.
  Eigen::Vector3f CenterOfLowerVoxel(const double x, const double y,
                                     const double z) const {
    // Center of the cell containing (x, y, z).
    Eigen::Vector3f center(static_cast<int>(x), static_cast<int>(y),
                                static_cast<int>(z));
    // Move to the next lower voxel center.
    if (center.x() > x) {
      center.x() -= 1;
    }
    if (center.y() > y) {
      center.y() -= 1;
    }
    if (center.z() > z) {
      center.z() -= 1;
    }
    return center;
  }

  const HybridDecayGrid& hybrid_grid_;
};

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_INTERPOLATED_DECAY_GRID_H_ */
