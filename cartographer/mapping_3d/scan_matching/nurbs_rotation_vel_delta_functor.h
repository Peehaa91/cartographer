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
//      Eigen::Quaterniond quat = transform::AngleAxisVectorToRotationQuaternion(rot_vel.second);
//      quat.normalize();
//      Eigen::Quaterniond quat_inv(quat.w(), -quat.x(), -quat.y(), -quat.z());

      initial_angular_vel_.push_back(std::make_pair(rot_vel.first,
                                                    rot_vel.second));
    }

  }

  NurbsRotationVelocityDeltaFunctor(const NurbsRotationVelocityDeltaFunctor&) = delete;
  NurbsRotationVelocityDeltaFunctor& operator=(const NurbsRotationVelocityDeltaFunctor&) = delete;

  template<typename T>
  double getScalar(const T& jet) const {
    return ceres::Jet<double, 35>(jet).a;
  }

  template<typename T>
  T getVelocityDiff(const T& v1, const double& v2) const
  {
    T diff = v1 - v2;
    if (diff < - T(M_PI)){
      diff += T(2*M_PI);
    }
    if (diff > T(M_PI)){
      diff -= T(2*M_PI);
    }
    return diff;
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

    int numSamplePoints = initial_angular_vel_.size();
      int counter = 0 ;
      for (int i = 0; i < numSamplePoints; ++i) {
        double point = ((common::ToSeconds(initial_angular_vel_[i].first - begin_)
                  / common::ToSeconds(end_ - begin_)));
        std::array<T, output_dim> output;
        nurbs_deriv.getPoint(&point, output.begin());
//        T delta[4];
//        T initial_rotation_inverse[4] = {
//            T(initial_angular_vel_[i].second.w()), T(initial_angular_vel_[i].second.x()),
//            T(initial_angular_vel_[i].second.y()), T(initial_angular_vel_[i].second.z())};
        Eigen::Quaternion<T> rotation_vel(
                    output[6], output[3],
                    output[4], output[5]);
        Eigen::Matrix<T, 3, 1> angular_vel = rotation_vel.toRotationMatrix().eulerAngles(0, 1, 2);
//        ceres::QuaternionProduct(initial_rotation_inverse, rotation_vel,
//                                 delta);

//        Eigen::Quaternion<T> quat = Eigen::Quaternion<T>(output[6],output[3],output[4],output[3]);
//        Eigen::Matrix<T, 3, 1> angles = quat.toRotationMatrix().eulerAngles(0, 1, 2);
        // Will compute the squared norm of the imaginary component of the delta
        // quaternion which is sin(phi/2)^2.
//        LOG(INFO)<<"output x:"<<getScalar(angular_vel[0])<<"init x: "<<initial_angular_vel_[i].second[0];
//        LOG(INFO)<<"output y:"<<getScalar(angular_vel[1])<<"init y: "<<initial_angular_vel_[i].second[1];
//        LOG(INFO)<<"output z:"<<getScalar(angular_vel[2])<<"init z: "<<initial_angular_vel_[i].second[2];
//        LOG(INFO)<<"diff x:"<<getScalar(getVelocityDiff(angular_vel[0], initial_angular_vel_[i].second[0]));
//        LOG(INFO)<<"diff y:"<<getScalar(getVelocityDiff(angular_vel[1], initial_angular_vel_[i].second[1]));
//        LOG(INFO)<<"diff z:"<<getScalar(getVelocityDiff(angular_vel[2], initial_angular_vel_[i].second[2]));
        residual[counter] = scaling_factor_ * getVelocityDiff(angular_vel[0], initial_angular_vel_[i].second[0]) ;
        counter++;
        residual[counter] = scaling_factor_ * getVelocityDiff(angular_vel[1], initial_angular_vel_[i].second[1]);
        counter++;
        residual[counter] = scaling_factor_ * getVelocityDiff(angular_vel[2], initial_angular_vel_[i].second[2]);
        counter++;
//        residual[counter] = scaling_factor_ * (angles[0] - T(initial_rotation_vel_inverse_[i].second[0]));
//        counter++;
//        residual[counter] = scaling_factor_ * (angles[1] - T(initial_rotation_vel_inverse_[i].second[1]));;
//        counter++;
//        residual[counter] = scaling_factor_ * (angles[2] - T(initial_rotation_vel_inverse_[i].second[2]));;
//        counter++;
      }
    return true;
  }

 private:
  const double scaling_factor_;
  std::vector<std::pair<common::Time, Eigen::Vector3d>> initial_angular_vel_;
  const common::Time begin_;
  const common::Time end_;
};

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_ROTATION_VEL_DELTA_FUNCTOR_H_ */
