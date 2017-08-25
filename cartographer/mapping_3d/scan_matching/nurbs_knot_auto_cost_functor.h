/*
 * nurbs_knot_auto_cost_functor.h
 *
 *  Created on: Aug 23, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_KNOT_AUTO_COST_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_KNOT_AUTO_COST_FUNCTOR_H_

#include "Eigen/Core"
#include "cartographer/mapping_3d/hybrid_grid.h"
#include "cartographer/mapping_3d/scan_matching/interpolated_grid.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/transform/transform.h"
#include "cartographer/mapping_3d/scan_matching/nurbs.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {
const int input_dim = 1;
const int output_dim = 7;
const int degree = 2;
const KnotType knot_type = KnotType::UNIFORM;
const WeightType weight_type = WeightType::RATIONAL;
// Computes the cost of inserting occupied space described by the point cloud
// into the map. The cost increases with the amount of free space that would be
// replaced by occupied space.
class NurbsKnotAutoCostFunctor {
 public:
  // Creates an OccupiedSpaceCostFunctor using the specified grid, 'rotation' to
  // add to all poses, and point cloud.
  NurbsKnotAutoCostFunctor(
      const double scaling_factor,
      const HybridGrid& hybrid_grid,
      const common::Time& begin,
      const common::Time& end,
      const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec)
      : scaling_factor_(scaling_factor),
        range_data_vec_(range_data_vec),
        interpolated_grid_(hybrid_grid),
        begin_(begin),
        end_(end) {
  }

  NurbsKnotAutoCostFunctor(const NurbsKnotCostFunctor&) = delete;
  NurbsKnotAutoCostFunctor& operator=(const NurbsKnotCostFunctor&) = delete;

  template <typename T>
  bool operator()(const T* const point_1_trans, const T* const point_1_rotation,
                  const T* const point_2_trans, const T* const point_2_rotation,
                  const T* const point_3_trans, const T* const point_3_rotation,
                  const T* const point_4_trans, const T* const point_4_rotation,
                  T* const residual) const {

    std::vector<const T*> vec_trans = {point_1_trans, point_2_trans,
        point_3_trans, point_4_trans};
    std::vector<const T*> vec_rot = {point_1_rotation, point_2_rotation,
        point_3_rotation, point_4_rotation};
    Nurbs<T, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<
          T, input_dim, output_dim, knot_type, weight_type>();
      NurbsKnots<T, input_dim, output_dim, knot_type> knots = NurbsKnots<
          T, input_dim, output_dim, knot_type>();
      NurbsWeights<T, input_dim, output_dim, weight_type> weights =
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
    points = nurbs.getPoints();
    int numSamplePoints = range_data_vec_.size();
    sensor::PointCloud complete_pointcloud;
    int counter = 0 ;
    for (int i = 0; i < numSamplePoints; ++i) {
  //    double point = minU
  //        + (maxU - minU) / (double) (numSamplePoints - 1) * (double) i;
      double point = ((common::ToSeconds(range_data_vec_[i].first - begin_)
          / common::ToSeconds(end_ - begin_)));
      std::array<T, output_dim> output;
      nurbs.getPoint(&point, output.begin());
      const transform::Rigid3<T> transform(
          Eigen::Matrix<T, 3, 1>(output[0], output[1], output[2]),
          Eigen::Quaternion<T>(output[6], output[3], output[4],
                               output[5]));
      for (Eigen::Vector3f point : range_data_vec_[i].second.returns)
      {
        Eigen::Matrix<T, 3, 1> world = transform* point.cast<T>();
        const T probability =
                    interpolated_grid_.GetProbability(world[0], world[1], world[2]);
                residual[counter] = scaling_factor_ * (1. - probability);
                counter++;
      }
    }
//    int number_of_residuals = 0;
//    for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec_)
//      number_of_residuals += data.second.returns.size();
//    for (int i = 0; i < number_of_residuals; i++)
//      residual[i] = T(0);
    return true;
  }
private:
  const double scaling_factor_;
  const InterpolatedGrid interpolated_grid_;
  const common::Time begin_;
  const common::Time end_;
  const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec_;
};


}
}
}




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_KNOT_AUTO_COST_FUNCTOR_H_ */
