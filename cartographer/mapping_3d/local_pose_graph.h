/*
 * local_pose_graph.h
 *
 *  Created on: Jul 15, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_LOCAL_POSE_GRAPH_H_
#define CARTOGRAPHER_MAPPING_3D_LOCAL_POSE_GRAPH_H_

#include "cartographer/mapping_3d/proto/local_pose_graph_options.pb.h"
#include "cartographer/mapping_3d/scan_matching/nurbs.h"
#include "cartographer/mapping/global_trajectory_builder_interface.h"
#include "cartographer/mapping_3d/local_trajectory_builder.h"
#include "cartographer/mapping_3d/scan_matching/ceres_scan_matcher.h"

namespace cartographer {
namespace mapping_3d {

const int input_dim = 1;
const int output_dim = 7;
const int degree = 3;
const KnotType knot_type = KnotType::UNIFORM;
const WeightType weight_type = WeightType::NON_RATIONAL;

proto::LocalPoseGraphOptions CreateLocalPoseGraphOptions(
    common::LuaParameterDictionary* parameter_dictionary);

struct PoseAndRangeData{

  PoseAndRangeData() = default;
  PoseAndRangeData(common::Time time, transform::Rigid3d pose, sensor::RangeData range_data) :
    time(time),
    pose(pose),
    range_data(range_data){}
  common::Time time = common::Time::min();
  transform::Rigid3d pose = transform::Rigid3d::Identity();
  sensor::RangeData range_data;
};
class LocalPoseGraph {
 public:
  using PoseEstimate = mapping::GlobalTrajectoryBuilderInterface::PoseEstimate;
  LocalPoseGraph(const proto::LocalPoseGraphOptions& options);
  std::vector<PoseAndRangeData> createSplineFromControlVector(
      std::vector<PoseEstimate>& pose_vec,
      std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
      std::vector<std::pair<common::Time, Eigen::Vector3d>>& imu_angular_vel_data_vector,
      std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector,
      const std::vector<const HybridGrid*>& hybrid_grid,
      const std::unique_ptr<scan_matching::CeresScanMatcher>& ceres_scan_matcher,
      const std::vector<const HybridDecayGrid*>& decay_grid);
  int nurbs_number = 0;

 private:

  const proto::LocalPoseGraphOptions options_;
  void createPoseAndScanFromSpline(std::vector<PoseEstimate>& control_point_vec,
                                 std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
                                 std::vector<PoseAndRangeData>& pose_and_cloud_vec);

  void createDerivativeSplines(std::vector<PoseEstimate>& control_point_vec,
                                std::vector<PoseAndRangeData> & range_data_vec,
                                std::vector<std::pair<common::Time, transform::Rigid3d>>& velocity_data,
                                std::vector<std::pair<common::Time, Eigen::Vector3d>>& angular_vel_data,
                                std::vector<std::pair<common::Time, transform::Rigid3d>>& acceleration_data,
                                std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector,
                                std::vector<std::pair<common::Time, Eigen::Vector3d>>& imu_res_vec);
  void writeSplineInFile(std::vector<PoseEstimate>& control_point_vec,
                         std::vector<PoseAndRangeData> & range_data_vec,
                         std::string file_name,
                         std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector,
                         std::vector<std::pair<common::Time, Eigen::Vector3d>>& imu_angular_vel_data_vector);
  /*B Spline Parameters*/
  unsigned int input_size_;
  unsigned int output_size_;
  KnotType knot_type_;
  WeightType weight_type_;
  Eigen::Vector3d last_lin_vel_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d last_ang_vel_ = Eigen::Vector3d::Zero();
//  Nurbs<double, input_size_, output_size_, knot_type_, weight_type_> nurbs_;

};
}
}

#endif /* CARTOGRAPHER_MAPPING_3D_LOCAL_POSE_GRAPH_H_ */
