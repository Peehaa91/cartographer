#ifndef CARTOGRAPHER_MAPPING_3D_CONTINUOUS_LOCAL_TRAJECTORY_BUILDER_H_
#define CARTOGRAPHER_MAPPING_3D_CONTINUOUS_LOCAL_TRAJECTORY_BUILDER_H_

#include "cartographer/common/time.h"
#include "cartographer/mapping/global_trajectory_builder_interface.h"
#include "cartographer/mapping/imu_tracker.h"
#include "cartographer/mapping/odometry_state_tracker.h"
#include "cartographer/mapping_3d/ground_plane_tracker.h"
#include "cartographer/mapping_3d/motion_filter.h"
#include "cartographer/mapping_3d/proto/local_trajectory_builder_options.pb.h"
#include "cartographer/mapping_3d/scan_matching/ceres_scan_matcher.h"
#include "cartographer/mapping_3d/scan_matching/real_time_correlative_scan_matcher.h"
#include "cartographer/mapping_3d/submaps.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/sensor/voxel_filter.h"
#include "cartographer/transform/rigid_transform.h"

namespace cartographer {
namespace mapping_3d {

class ContinuousLocalTrajectoryBuilder {
 public:
  using PoseEstimate = mapping::GlobalTrajectoryBuilderInterface::PoseEstimate;
  ContinuousLocalTrajectoryBuilder();

  explicit ContinuousLocalTrajectoryBuilder(
      const proto::LocalTrajectoryBuilderOptions& options);

  void AddImuData(common::Time time, const Eigen::Vector3d& linear_acceleration,
                  const Eigen::Vector3d& angular_velocity);
  void AddRangefinderData(common::Time time, const Eigen::Vector3f& origin,
                          const sensor::PointCloud& ranges);
  void AddOdometerData(common::Time time,
                       const transform::Rigid3d& odometer_pose);

  void AddPlaneData(const common::Time time,
                    const Eigen::Vector4d& coefficients);

 private:
  std::vector<PoseEstimate> pose_vec_;
};
}
}

#endif
