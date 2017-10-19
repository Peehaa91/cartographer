/*
 * nurbs_translation_acceleration_delta_functor.h
 *
 *  Created on: Sep 28, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_TRANSLATION_ACCELERATION_DELTA_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_TRANSLATION_ACCELERATION_DELTA_FUNCTOR_H_

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
class NurbsTranslationAccelerationDeltaFunctor {
 public:
  // Constructs a new RotationDeltaCostFunctor from the given
  // 'rotation_estimate'.
  explicit NurbsTranslationAccelerationDeltaFunctor(const double scaling_factor,
                                    const std::vector<std::pair<common::Time,
                                    Eigen::Vector3d>>& initial_linear_acceleration,
                                    const common::Time& begin,
                                    const common::Time& end)
      : scaling_factor_(scaling_factor),
        initial_linear_acceleration_(initial_linear_acceleration),
        begin_(begin),
        end_(end){
  }

  NurbsTranslationAccelerationDeltaFunctor(const NurbsTranslationAccelerationDeltaFunctor&) = delete;
  NurbsTranslationAccelerationDeltaFunctor& operator=(const NurbsTranslationAccelerationDeltaFunctor&) = delete;

  template<typename T>
  double getScalar(const T& jet) const {
    return ceres::Jet<double, 35>(jet).a;
  }

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
    Nurbs<T, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<
          T, input_dim, output_dim, knot_type, weight_type>();
    Nurbs<T, input_dim, output_dim, knot_type, weight_type> nurbs_second_deriv = Nurbs<
          T, input_dim, output_dim, knot_type, weight_type>();
    NurbsKnots<T, input_dim, output_dim, knot_type> knots = NurbsKnots<
        T, input_dim, output_dim, knot_type>();
    NurbsKnots<T, input_dim, output_dim, knot_type> knots_first_deriv = NurbsKnots<
        T, input_dim, output_dim, knot_type>();
    NurbsKnots<T, input_dim, output_dim, knot_type> knots_second_deriv = NurbsKnots<
        T, input_dim, output_dim, knot_type>();
    NurbsWeights<T, input_dim, output_dim, weight_type> weights =
        NurbsWeights<T, input_dim, output_dim, weight_type>();
    NurbsWeights<T, input_dim, output_dim, weight_type> weights_second_deriv =
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
    weights.create(points.shape());
    nurbs.init(degree, knots, weights, points);

    std::array<std::vector<T>, input_dim> knot_vec = knots.getKnot();

    boost::multi_array<std::array<T, output_dim>, input_dim> points_first_deriv(
        boost::extents[vec_trans.size() -1]);
    for (int i = 0; i < vec_trans.size() -1; ++i) {
      points_first_deriv[i][0] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][0] - points[i][0]);
      points_first_deriv[i][1] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][1] - points[i][1]);
      points_first_deriv[i][2] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][2] - points[i][2]);
      points_first_deriv[i][3] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][3] - points[i][3]);
      points_first_deriv[i][4] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][4] - points[i][4]);
      points_first_deriv[i][5] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][5] - points[i][5]);
      points_first_deriv[i][6] = T(degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][6] - points[i][6]);
    }
    knots_first_deriv.create(degree - 1, points_first_deriv.shape());



    std::array<std::vector<T>, input_dim> knot_first_vec = knots_first_deriv.getKnot();

    boost::multi_array<std::array<T, output_dim>, input_dim> points_second_deriv(
        boost::extents[vec_trans.size() -2]);
    for (int i = 0; i < vec_trans.size() -2; ++i) {
      points_second_deriv[i][0] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][0] - points_first_deriv[i][0]);
      points_second_deriv[i][1] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][1] - points_first_deriv[i][1]);
      points_second_deriv[i][2] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][2] - points_first_deriv[i][2]);
      points_second_deriv[i][3] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][3] - points_first_deriv[i][3]);
      points_second_deriv[i][4] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][4] - points_first_deriv[i][4]);
      points_second_deriv[i][5] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][5] - points_first_deriv[i][5]);
      points_second_deriv[i][6] = T(degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
          * (points_first_deriv[i + 1][6] - points_first_deriv[i][6]);
    }

    weights_second_deriv.create(points_second_deriv.shape());
    knots_second_deriv.create(degree - 2, points_second_deriv.shape());
    nurbs_second_deriv.init(degree - 2, knots_second_deriv, weights_second_deriv, points_second_deriv);
    int numSamplePoints = initial_linear_acceleration_.size();
    int counter = 0 ;
    for (int i = 0; i < numSamplePoints; ++i) {
      double point = ((common::ToSeconds(initial_linear_acceleration_[i].first - begin_)
                / common::ToSeconds(end_ - begin_)));
      std::array<T, output_dim> output_second_deriv;
      std::array<T, output_dim> output;
      nurbs_second_deriv.getPoint(&point, output_second_deriv.begin());
      nurbs.getPoint(&point, output.begin());
      Eigen::Quaternion<T>rotation(output[6], output[3], output[4], output[5]);
      Eigen::Matrix<T, 3, 1> lin_acc(output_second_deriv[0], output_second_deriv[1], output_second_deriv[2]);
      Eigen::Matrix<T, 3, 1> gravity(T(0), T(0), T(9.81));
      Eigen::Matrix<T, 3, 1> gravity_rotated = rotation.inverse() * gravity;
      Eigen::Matrix<T, 3, 1> inital_acceleration = initial_linear_acceleration_[i].second.cast<T>() - gravity_rotated;

//      LOG(INFO)<<"output x:"<<getScalar(lin_acc[0])<<"init x: "<<getScalar(inital_acceleration[0]);
//      LOG(INFO)<<"output y:"<<getScalar(lin_acc[1])<<"init y: "<<getScalar(inital_acceleration[1]);
//      LOG(INFO)<<"output z:"<<getScalar(lin_acc[2])<<"init z: "<<getScalar(inital_acceleration[2]);

      residual[counter] = scaling_factor_ * (inital_acceleration[0] - lin_acc[0]);
      counter++;
      residual[counter] = scaling_factor_ * (inital_acceleration[1] - lin_acc[1]);
      counter++;
      residual[counter] = scaling_factor_ * (inital_acceleration[2] - lin_acc[2]);
      counter++;
    }
    return true;
  }

 private:
  const double scaling_factor_;
  std::vector<std::pair<common::Time, Eigen::Vector3d>> initial_linear_acceleration_;
  const common::Time begin_;
  const common::Time end_;
};

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_TRANSLATION_ACCELERATION_DELTA_FUNCTOR_H_ */
