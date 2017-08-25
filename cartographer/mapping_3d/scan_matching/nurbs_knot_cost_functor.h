/*
 * nurbs_knot_cost_functor.h
 *
 *  Created on: Aug 14, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_KNOT_COST_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_KNOT_COST_FUNCTOR_H_

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
class NurbsKnotCostFunctor {
 public:
  // Creates an OccupiedSpaceCostFunctor using the specified grid, 'rotation' to
  // add to all poses, and point cloud.
  NurbsKnotCostFunctor(
      const double scaling_factor,
      const HybridGrid& hybrid_grid,
      const common::Time& begin,
      const common::Time& end,
      const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec,
      Nurbs<double, 1, 7, KnotType::UNIFORM, WeightType::RATIONAL>* nurbs)
      : scaling_factor_(scaling_factor),
        range_data_vec_(range_data_vec),
        interpolated_grid_(hybrid_grid),
        begin_(begin),
        end_(end),
        nurbs_(nurbs) {
  }

  NurbsKnotCostFunctor(const NurbsKnotCostFunctor&) = delete;
  NurbsKnotCostFunctor& operator=(const NurbsKnotCostFunctor&) = delete;

  bool operator()(const double* const point_1, const double* const point_2,
                  const double* const point_3, const double* const point_4,
                  const double* const point_5, double* residual) const {
//    LOG(INFO)<<"point_1: "<<point_1[0]<<","<<point_1[1]<<","<<point_1[2]<<","<<point_1[3]<<","<<point_1[4]<<","<<point_1[5];
//    LOG(INFO)<<"point_2: "<<point_2[0]<<","<<point_2[1]<<","<<point_2[2]<<","<<point_2[3]<<","<<point_2[4]<<","<<point_2[5];
//    LOG(INFO)<<"point_3: "<<point_3[0]<<","<<point_3[1]<<","<<point_3[2]<<","<<point_3[3]<<","<<point_3[4]<<","<<point_3[5];
//    LOG(INFO)<<"point_4: "<<point_4[0]<<","<<point_4[1]<<","<<point_4[2]<<","<<point_4[3]<<","<<point_4[4]<<","<<point_4[5];
//    LOG(INFO)<<"point_5: "<<point_5[0]<<","<<point_5[1]<<","<<point_5[2]<<","<<point_5[3]<<","<<point_5[4]<<","<<point_5[5];

    std::vector<const double*> vec = {point_1, point_2, point_3, point_4, point_5};
    boost::multi_array<std::array<double, 7>, 1>& points = nurbs_->getPoints();
    for (int i = 0; i < vec.size(); i++) {
      points[i][0] = vec[i][0];
      points[i][1] = vec[i][1];
      points[i][2] = vec[i][2];
      points[i][3] = vec[i][3];
      points[i][4] = vec[i][4];
      points[i][5] = vec[i][5];
      points[i][6] = vec[i][6];

    }
    int numSamplePoints = range_data_vec_.size();
    sensor::PointCloud complete_pointcloud;
    int counter = 0 ;
    for (int i = 0; i < numSamplePoints; ++i) {
  //    double point = minU
  //        + (maxU - minU) / (double) (numSamplePoints - 1) * (double) i;
      double point = common::ToSeconds(range_data_vec_[i].first - begin_)
          / common::ToSeconds(end_ - begin_);
      std::array<double, 7> output;
      nurbs_->getPoint(&point, output.begin());
      transform::Rigid3d pose = transform::Rigid3d(
          Eigen::Vector3d(output[0], output[1], output[2]),
          Eigen::Quaterniond(output[6], output[3], output[4], output[5]));
//          transform::AngleAxisVectorToRotationQuaternion(
//              Eigen::Vector3d(output[3], output[4], output[5])));
//      sensor::PointCloud pointcloud = sensor::TransformPointCloud(range_data_vec_[i].second.returns,
      for (Eigen::Vector3f point : range_data_vec_[i].second.returns)
      {
        Eigen::Matrix<double, 3, 1> world = pose* point.cast<double>();
        const float probability =
                    interpolated_grid_.GetProbability(world[0], world[1], world[2]);
                residual[counter] = scaling_factor_ * (1. - probability);
                counter++;
      }
//                                        pose.cast<float>());
//      for (size_t j = 0; j < pointcloud.size(); ++j) {
//        const Eigen::Matrix<double, 3, 1> world =
//            pointcloud[j].cast<double>();
//        const double probability =
//            interpolated_grid_.GetProbability(world[0], world[1], world[2]);
//        residual[counter] = scaling_factor_ * (1. - probability);
//        counter++;
//      }
    }
//    int number_of_residuals = 0;
//    for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec_)
//      number_of_residuals += data.second.returns.size();
//    for (int i = 0; i < number_of_residuals; i++)
//      residual[i] = 0;
    return true;
  }
private:
  const double scaling_factor_;
  const InterpolatedGrid interpolated_grid_;
  const common::Time begin_;
  const common::Time end_;
  const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec_;
  Nurbs<double, 1, 7, KnotType::UNIFORM, WeightType::RATIONAL>* nurbs_;
};
}
}
}

#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_KNOT_COST_FUNCTOR_H_ */
