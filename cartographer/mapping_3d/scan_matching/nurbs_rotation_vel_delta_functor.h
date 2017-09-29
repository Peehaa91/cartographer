/*
 * nurbs_rotation_vel_delta_functor.h
 *
 *  Created on: Sep 28, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_ROTATION_VEL_DELTA_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_ROTATION_VEL_DELTA_FUNCTOR_H_

#include <cmath>

#include "Eigen/Core"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/transform/transform.h"
#include "ceres/rotation.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_control_point_auto_cost_functor.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {

// Computes the cost of rotating the pose estimate. Cost increases with the
// solution's distance from the rotation estimate.
class NurbsRotationVelocityDeltaFunctor {
 public:
  // Constructs a new RotationDeltaCostFunctor from the given
  // 'rotation_estimate'.
  explicit NurbsRotationVelocityDeltaFunctor(const double scaling_factor,
                                    const std::vector<std::pair<common::Time,
                                    Eigen::Vector3d>>& initial_rotation,
                                    const common::Time& begin,
                                    const common::Time& end)
      : scaling_factor_(scaling_factor),
        begin_(begin),
        end_(end){
    for (auto& rot_vel : initial_rotation)
    {
      Eigen::Quaterniond quat = transform::AngleAxisVectorToRotationQuaternion(rot_vel.second);
      Eigen::Quaterniond quat_inv(quat.w(), -quat.x(), -quat.y(), -quat.z());
      initial_rotation_vel_inverse_.push_back(std::make_pair(rot_vel.first,
                                                             quat_inv));
    }

  }

  NurbsRotationVelocityDeltaFunctor(const NurbsRotationVelocityDeltaFunctor&) = delete;
  NurbsRotationVelocityDeltaFunctor& operator=(const NurbsRotationVelocityDeltaFunctor&) = delete;

  template <typename T>
  bool operator()(const T* const point_1_trans, const T* const point_1_rotation,
                  const T* const point_2_trans, const T* const point_2_rotation,
                  const T* const point_3_trans, const T* const point_3_rotation,
                  const T* const point_4_trans, const T* const point_4_rotation,
                  const T* const point_5_trans, const T* const point_5_rotation,
                  T* const residual) const {

    std::vector<const T*> vec_trans = {point_1_trans, point_2_trans,
            point_3_trans, point_4_trans, point_5_trans};
    std::vector<const T*> vec_rot = {point_1_rotation, point_2_rotation,
        point_3_rotation, point_4_rotation, point_5_rotation};
    Nurbs<T, input_dim, output_dim, knot_type, weight_type> nurbs_deriv = Nurbs<
          T, input_dim, output_dim, knot_type, weight_type>();
    NurbsKnots<T, input_dim, output_dim, knot_type> knots = NurbsKnots<
        T, input_dim, output_dim, knot_type>();
    NurbsKnots<T, input_dim, output_dim, knot_type> knots_deriv = NurbsKnots<
        T, input_dim, output_dim, knot_type>();
    NurbsWeights<T, input_dim, output_dim, weight_type> weights_deriv =
        NurbsWeights<T, input_dim, output_dim, weight_type>();
    boost::multi_array<std::array<T, output_dim>, input_dim> points( boost::extents[vec_trans.size()]);
    for (int i = 0; i < vec_trans.size(); i++) {
      points[i][0] = vec_trans[i][0];
      points[i][1] = vec_trans[i][1];
      points[i][2] = vec_trans[i][2];
      points[i][3] = vec_rot[i][1];
      points[i][4] = vec_rot[i][2];
      points[i][5] = vec_rot[i][3];
      points[i][6] = vec_rot[i][0];
    }
    knots.create(degree, points.shape());

    std::array<std::vector<T>, input_dim> knot_vec = knots.getKnot();

    boost::multi_array<std::array<T, output_dim>, input_dim> points_deriv(
        boost::extents[vec_trans.size() -1]);
    for (int i = 0; i < vec_trans.size() -1; ++i) {
      points_deriv[i][0] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][0] - points[i][0]);
      points_deriv[i][1] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][1] - points[i][1]);
      points_deriv[i][2] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][2] - points[i][2]);
      points_deriv[i][3] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][3] - points[i][3]);
      points_deriv[i][4] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][4] - points[i][4]);
      points_deriv[i][5] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][5] - points[i][5]);
      points_deriv[i][6] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][6] - points[i][6]);
    }
    knots_deriv.create(degree - 1, points_deriv.shape());
    weights_deriv.create(points_deriv.shape());
    nurbs_deriv.init(degree -1, knots_deriv, weights_deriv, points_deriv);

    int numSamplePoints = initial_rotation_vel_inverse_.size();
      int counter = 0 ;
      for (int i = 0; i < numSamplePoints; ++i) {
        double point = ((common::ToSeconds(initial_rotation_vel_inverse_[i].first - begin_)
                  / common::ToSeconds(end_ - begin_)));
        std::array<T, output_dim> output;
        nurbs_deriv.getPoint(&point, output.begin());
        T delta[4];
        T initial_rotation_inverse[4] = {
            T(initial_rotation_vel_inverse_[i].second.w()), T(initial_rotation_vel_inverse_[i].second.x()),
            T(initial_rotation_vel_inverse_[i].second.y()), T(initial_rotation_vel_inverse_[i].second.z())};
        T rotation_vel[4] = {
                    T(output[6]), T(output[3]),
                    T(output[4]), T(output[5])};
        ceres::QuaternionProduct(initial_rotation_inverse, rotation_vel,
                                 delta);
        // Will compute the squared norm of the imaginary component of the delta
        // quaternion which is sin(phi/2)^2.
        residual[counter] = scaling_factor_ * delta[1];
        counter++;
        residual[counter] = scaling_factor_ * delta[2];
        counter++;
        residual[counter] = scaling_factor_ * delta[3];
        counter++;
      }
    return true;
  }

 private:
  const double scaling_factor_;
  std::vector<std::pair<common::Time, Eigen::Quaterniond>> initial_rotation_vel_inverse_;
  const common::Time begin_;
  const common::Time end_;
};

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_ROTATION_VEL_DELTA_FUNCTOR_H_ */
