/*
 * nurbs_weight_cost_functor.h
 *
 *  Created on: Aug 8, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_WEIGHT_COST_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_WEIGHT_COST_FUNCTOR_H_

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
class NurbsWeightCostFunctor {
 public:
  // Creates an OccupiedSpaceCostFunctor using the specified grid, 'rotation' to
  // add to all poses, and point cloud.
  NurbsWeightCostFunctor(const double scaling_factor,
                           const HybridGrid& hybrid_grid,
                           const common::Time& begin,
                           const common::Time& end,
                           const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec,
                           Nurbs<double, 1, 6, KnotType::UNIFORM, WeightType::RATIONAL>* nurbs)
      : scaling_factor_(scaling_factor),
        range_data_vec_(range_data_vec),
        interpolated_grid_(hybrid_grid),
        begin_(begin),
        end_(end),
        nurbs_(nurbs){}

  NurbsWeightCostFunctor(const NurbsWeightCostFunctor&) = delete;
  NurbsWeightCostFunctor& operator=(const NurbsWeightCostFunctor&) = delete;

//  template <typename T>
  bool operator()(const double* const weights,
                  double* residual) const {
//    LOG(INFO)<<"optimization called";
//    Nurbs<double, 1, 6, KnotType::UNIFORM, WeightType::RATIONAL> nurb_copy(nurbs_);
    NurbsWeights<double, 1, 6, WeightType::RATIONAL>& weight = nurbs_->getWeights();
    boost::multi_array<double, 1>& weights_vec =  weight.getWeights();
    for (int i = 0; i < 5; i++)
    {
      weights_vec[i] = weights[i];
//      LOG(INFO)<<"weight in opt: "<<weights[i];
    }
//    double minU = nurbs_.getMinU(0);
//    double maxU = nurbs_.getMaxU(0);

    int numSamplePoints = range_data_vec_.size();
    sensor::PointCloud complete_pointcloud;
    int counter = 0 ;
    for (int i = 0; i < numSamplePoints; ++i) {
  //    double point = minU
  //        + (maxU - minU) / (double) (numSamplePoints - 1) * (double) i;
      double point = common::ToSeconds(range_data_vec_[i].first - begin_)
          / common::ToSeconds(end_ - begin_);
      std::array<double, 6> output;
      nurbs_->getPoint(&point, output.begin());
      transform::Rigid3d pose = transform::Rigid3d(
          Eigen::Vector3d(output[0], output[1], output[2]),
          transform::AngleAxisVectorToRotationQuaternion(
              Eigen::Vector3d(output[3], output[4], output[5])));
      sensor::PointCloud pointcloud = sensor::TransformPointCloud(range_data_vec_[i].second.returns,
                                        pose.cast<float>());
      for (size_t j = 0; j < pointcloud.size(); ++j) {
        const Eigen::Matrix<double, 3, 1> world =
            pointcloud[j].cast<double>();
        const double probability =
            interpolated_grid_.GetProbability(world[0], world[1], world[2]);
        residual[counter] = scaling_factor_ * (1. - probability);
        counter++;
      }
//      complete_pointcloud.insert(complete_pointcloud.end(), pointcloud.begin(), pointcloud.end());
    }

//    for (size_t i = 0; i < complete_pointcloud.size(); ++i) {
//      const Eigen::Matrix<double, 3, 1> world =
//          complete_pointcloud[i].cast<double>();
//      const double probability =
//          interpolated_grid_.GetProbability(world[0], world[1], world[2]);
//      residual[i] = scaling_factor_ * (1. - probability);
//    }
    return true;
  }

 private:
  const double scaling_factor_;
  const InterpolatedGrid interpolated_grid_;
  const common::Time begin_;
  const common::Time end_;
  const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec_;
  Nurbs<double, 1, 6, KnotType::UNIFORM, WeightType::RATIONAL>* nurbs_;
};

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_WEIGHT_COST_FUNCTOR_H_ */
