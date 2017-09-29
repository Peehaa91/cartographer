#ifndef CARTOGRAPHER_MAPPING_3D_CONTINUOUS_LOCAL_TRAJECTORY_BUILDER_H_
#define CARTOGRAPHER_MAPPING_3D_CONTINUOUS_LOCAL_TRAJECTORY_BUILDER_H_

#include "cartographer/common/time.h"
#include "cartographer/mapping/global_trajectory_builder_interface.h"
#include "cartographer/mapping/imu_tracker.h"
#include "cartographer/mapping/odometry_state_tracker.h"
#include "cartographer/mapping_3d/ground_plane_tracker.h"
#include "cartographer/mapping_3d/motion_filter.h"
#include "cartographer/mapping_3d/proto/local_trajectory_builder_options.pb.h"
#include "cartographer/mapping_3d/local_trajectory_builder.h"
#include "cartographer/mapping_3d/scan_matching/ceres_scan_matcher.h"
#include "cartographer/mapping_3d/scan_matching/real_time_correlative_scan_matcher.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/mapping_3d/local_pose_graph.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/sensor/voxel_filter.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/sensor/voxel_filter.h"

namespace cartographer {
namespace mapping_3d {

class ContinuousLocalTrajectoryBuilder {
 public:
  using PoseEstimate = mapping::GlobalTrajectoryBuilderInterface::PoseEstimate;

  explicit ContinuousLocalTrajectoryBuilder(
      const proto::LocalTrajectoryBuilderOptions& options);

  void AddImuData(common::Time time, const Eigen::Vector3d& linear_acceleration,
                  const Eigen::Vector3d& angular_velocity);
  std::unique_ptr<LocalTrajectoryBuilder::InsertionResult> AddRangefinderData(common::Time time, const Eigen::Vector3f& origin,
                          const sensor::PointCloud& ranges);
  void AddOdometerData(common::Time time,
                       const transform::Rigid3d& odometer_pose);

  void AddPlaneData(const common::Time time,
                    const Eigen::Vector4d& coefficients);

  const mapping::GlobalTrajectoryBuilderInterface::PoseEstimate& pose_estimate() const;

 private:

  std::unique_ptr<LocalTrajectoryBuilder::InsertionResult> InsertIntoSubmap(
       common::Time time, const sensor::RangeData& range_data_in_tracking,
       const transform::Rigid3d& pose_observation);

  std::unique_ptr<LocalTrajectoryBuilder::InsertionResult> AddAccumulatedRangeData(
        common::Time time, const sensor::RangeData& range_data_in_tracking);

  void Predict(common::Time time);
  const proto::LocalTrajectoryBuilderOptions options_;
  bool initial_imu_ = true;
  bool first_scan_ = true;
  common::Time last_control_point_time_;
  common::Time last_update_time_ = common::Time::min();
  common::Time last_scan_match_time_ = common::Time::min();
  transform::Rigid3d last_control_pose_estimate_ = transform::Rigid3d::Identity();
  mapping_3d::ActiveSubmapsDecay active_submaps_;
  std::unique_ptr<scan_matching::RealTimeCorrelativeScanMatcher>
      real_time_correlative_scan_matcher_;

  std::unique_ptr<scan_matching::CeresScanMatcher> ceres_scan_matcher_;
  Eigen::Vector3d velocity_estimate_ = Eigen::Vector3d::Zero();
  mapping::GlobalTrajectoryBuilderInterface::PoseEstimate last_pose_estimate_;// = {common::Time(),transform::Rigid3d::Identity(), {}};
  std::vector<PoseEstimate> pose_vec_;
  std::unique_ptr<mapping::ImuTracker> imu_tracker_;
  std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vector_;
  std::vector<std::pair<common::Time, Eigen::Vector3d>> imu_angular_vel_data_vector_;
  std::vector<std::pair<common::Time, Eigen::Vector3d>> linear_acc_data_vector_;
  Eigen::Vector3d acceleration_estimate_ = Eigen::Vector3d::Zero();
  LocalPoseGraph local_pose_graph_;
  int num_accumulated_ = 0;
  int control_point_counter_ = 0;
  Eigen::Vector3d vel_est_;
  transform::Rigid3f first_pose_estimate_ = transform::Rigid3f::Identity();
  transform::Rigid3d pose_estimate_ = transform::Rigid3d::Identity();
  sensor::RangeData accumulated_range_data_;
  std::shared_ptr<mapping_3d::SubmapDecay> scan_matcher_submap_;

};
}
}

#endif
