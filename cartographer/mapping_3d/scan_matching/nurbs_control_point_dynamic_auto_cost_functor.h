/*
 * nurbs_control_point_dynamic_auto_functor.h
 *
 *  Created on: Sep 7, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_CONTROL_POINT_DYNAMIC_AUTO_COST_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_CONTROL_POINT_DYNAMIC_AUTO_COST_FUNCTOR_H_

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
// Computes the cost of inserting occupied space described by the point cloud
// into the map. The cost increases with the amount of free space that would be
// replaced by occupied space.
class NurbsControlPointDynamicCostFunctor {
 public:
  // Creates an OccupiedSpaceCostFunctor using the specified grid, 'rotation' to
  // add to all poses, and point cloud.
  NurbsControlPointDynamicCostFunctor(
      const double scaling_factor,
      const int number_of_control_points,
      const HybridGrid& hybrid_grid,
      const common::Time& begin,
      const common::Time& end,
      const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec)
      : scaling_factor_(scaling_factor),
        number_of_control_points_(number_of_control_points),
        range_data_vec_(range_data_vec),
        interpolated_grid_(hybrid_grid),
        begin_(begin),
        end_(end) {
  }

  NurbsControlPointDynamicCostFunctor(const NurbsControlPointDynamicCostFunctor&) = delete;
  NurbsControlPointDynamicCostFunctor& operator=(const NurbsControlPointDynamicCostFunctor&) = delete;

  template <typename T>
  bool operator()(T const* const* parameters, T* residuals) const {

    std::vector<T> vec_trans;
    std::vector<T> vec_rot;
    Nurbs<T, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<
          T, input_dim, output_dim, knot_type, weight_type>();
      NurbsKnots<T, input_dim, output_dim, knot_type> knots = NurbsKnots<
          T, input_dim, output_dim, knot_type>();
      NurbsWeights<T, input_dim, output_dim, weight_type> weights =
          NurbsWeights<T, input_dim, output_dim, weight_type>();
    boost::multi_array<std::array<T, output_dim>, input_dim> points( boost::extents[number_of_control_points_]);
    int j = 0;
    for (int i = 0; i < 7*number_of_control_points_; i++){
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][0] = (*parameters)[i];
      i++;
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][1] = (*parameters)[i];
      i++;
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][2] = (*parameters)[i];
      i++;
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][6] = (*parameters)[i];
      i++;
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][3] = (*parameters)[i];
      i++;
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][4] = (*parameters)[i];
      i++;
      LOG(INFO)<<"i: "<<i<<" value: "<<(*parameters)[i];
      points[j][5] = (*parameters)[i];
      j++;
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
                    residuals[counter] = scaling_factor_ * (1. - probability);
                    counter++;
          }
        }
//    int number_of_residuals = 0;
//    for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec_)
//      number_of_residuals += data.second.returns.size();
//    for (int i = 0; i < number_of_residuals; i++)
//      residuals[i] = T(0);
    return true;
  }
private:
  const double scaling_factor_;
  const int number_of_control_points_;
  const InterpolatedGrid interpolated_grid_;
  const common::Time begin_;
  const common::Time end_;
  const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec_;
};


}
}
}



#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_CONTROL_POINT_DYNAMIC_AUTO_COST_FUNCTOR_H_ */
