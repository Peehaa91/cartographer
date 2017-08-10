#include <cartographer/mapping_3d/continuous_local_trajectory_builder.h>

namespace cartographer {
namespace mapping_3d {

ContinuousLocalTrajectoryBuilder::ContinuousLocalTrajectoryBuilder(
    const proto::LocalTrajectoryBuilderOptions& options)
    : options_(options),
      active_submaps_(options.submaps_options()),
      real_time_correlative_scan_matcher_(
          common::make_unique<scan_matching::RealTimeCorrelativeScanMatcher>(
              options_.real_time_correlative_scan_matcher_options())),
      ceres_scan_matcher_(
          common::make_unique<scan_matching::CeresScanMatcher>(
              options_.ceres_scan_matcher_options())),
      local_pose_graph_(LocalPoseGraph(options.local_pose_graph_options())) {
  LOG(INFO)<<"continuous builder created";
}

double counter = 1;
void ContinuousLocalTrajectoryBuilder::AddImuData(
    common::Time time, const Eigen::Vector3d& linear_acceleration,
    const Eigen::Vector3d& angular_velocity) {
  if (initial_imu_) {
    imu_tracker_ = common::make_unique<mapping::ImuTracker>(
        options_.imu_gravity_time_constant(), time);
  }
  const Eigen::Quaterniond last_orientation = imu_tracker_->orientation();
  //LOG(INFO)<<"gravity vec: "<<imu_tracker_->gravity_vector();
//  imu_tracker_->Advance(time);
  Predict(time);
  imu_tracker_->AddImuLinearAccelerationObservation(linear_acceleration);
  imu_tracker_->AddImuAngularVelocityObservation(angular_velocity);
  if (initial_imu_) {
    LOG(INFO)<<"initial imu";
//    last_update_time_ = time;
    last_control_point_time_ = time;
    pose_estimate_ = {transform::Rigid3d::Rotation(imu_tracker_->orientation())};
    acceleration_estimate_ = imu_tracker_->linear_acceleration();
    initial_imu_= false;
    return;
  }
//  if ((last_control_point_time_ - time) < common::Time::duration(1))
//  {
//    pose_vec_.push_back(imu_tracker_->)
//  }
//  const double delta_t = common::ToSeconds(time - last_update_time_);
//  acceleration_estimate_mean_ = (counter - 1)/counter * acceleration_estimate_mean_ + 1/counter * imu_tracker_->linear_acceleration();
//  counter++;
////  Eigen::Vector3d acceleration_estimate = (imu_tracker_->orientation()*linear_acceleration - imu_tracker_->gravity_vector());
////  LOG(INFO)<<"acceleration "<<acceleration_estimate<<std::endl<<" linear acceleration: "<<linear_acceleration<<std::endl<<" gravity vector: "<<imu_tracker_->gravity_vector();
//  // Constant velocity model.
//  const Eigen::Vector3d translation = last_pose_estimate_.pose.translation()
//      + delta_t * velocity_estimate_;// + 0.5 * acceleration_estimate * delta_t * delta_t;
//  velocity_estimate_ += acceleration_estimate_ * delta_t;
//  const Eigen::Quaterniond rotation = last_pose_estimate_.pose.rotation() *
//      last_orientation.inverse() *
//      imu_tracker_->orientation();
//  transform::Rigid3d pose_estimate = transform::Rigid3d(translation, rotation);
//  last_pose_estimate_ = {time, pose_estimate,{}};
//  last_update_time_ = time;
//  LOG(INFO)<<"pose: "<<last_pose_estimate_.pose<<" vel: "<<velocity_estimate_;

}

void ContinuousLocalTrajectoryBuilder::AddOdometerData(
    common::Time time, const transform::Rigid3d& odometer_pose) {

}

void ContinuousLocalTrajectoryBuilder::AddPlaneData(
    const common::Time time, const Eigen::Vector4d& coefficients) {

}

void ContinuousLocalTrajectoryBuilder::Predict(common::Time time) {
  CHECK(imu_tracker_ != nullptr);
  const Eigen::Quaterniond last_orientation = imu_tracker_->orientation();
  imu_tracker_->Advance(time);
  if (last_update_time_ > common::Time::min()) {
    const double delta_t = common::ToSeconds(time - last_update_time_);
    // Constant velocity model.
    const Eigen::Vector3d translation = pose_estimate_.translation()
        + delta_t * velocity_estimate_;
    const Eigen::Quaterniond rotation = pose_estimate_.rotation()
        * last_orientation.inverse() * imu_tracker_->orientation();
    transform::Rigid3d pose = transform::Rigid3d(translation, rotation);
//    acceleration_estimate_ = imu_tracker_->linear_acceleration();
//    velocity_estimate_ += acceleration_estimate_ * delta_t;
//      last_pose_estimate_ = PoseEstimate(time, pose,  {});
    pose_estimate_ = pose;
  }
  last_update_time_ = time;
}

std::unique_ptr<LocalTrajectoryBuilder::InsertionResult> ContinuousLocalTrajectoryBuilder::AddRangefinderData(
    const common::Time time, const Eigen::Vector3f& origin,
    const sensor::PointCloud& ranges) {
  if (imu_tracker_ == nullptr) {
    LOG(INFO)<< "ImuTracker not yet initialized.";
    return nullptr;
  }
  Predict(time);
  if ((range_data_vector_.empty() && pose_vec_.empty()) || control_point_counter_ == 32){
    pose_vec_.push_back( {time, pose_estimate_, {}});
    control_point_counter_ = 0;
  }
  sensor::AdaptiveVoxelFilter adaptive_voxel_filter(
      options_.high_resolution_adaptive_voxel_filter_options());
  const sensor::PointCloud filtered_point_cloud_in_tracking =
  adaptive_voxel_filter.Filter(ranges);
  control_point_counter_++;
  std::pair<common::Time, sensor::RangeData> range_data = std::make_pair(time, sensor::RangeData {origin, filtered_point_cloud_in_tracking, {}});
  range_data_vector_.push_back(range_data);
  if (num_accumulated_ == 0) {
    first_pose_estimate_ = pose_estimate_.cast<float>();
    accumulated_range_data_ =
    sensor::RangeData {Eigen::Vector3f::Zero(), {}, {}};
  }

  const transform::Rigid3f tracking_delta =
  first_pose_estimate_.inverse() * pose_estimate_.cast<float>();
  const sensor::RangeData range_data_in_first_tracking =
  sensor::TransformRangeData(sensor::RangeData {origin, ranges, {}},
      tracking_delta);
  for (const Eigen::Vector3f& hit : range_data_in_first_tracking.returns) {
    const Eigen::Vector3f delta = hit - range_data_in_first_tracking.origin;
    const float range = delta.norm();
    if (range >= options_.min_range()) {
      if (range <= options_.max_range()) {
        accumulated_range_data_.returns.push_back(hit);
      } else {
        // We insert a ray cropped to 'max_range' as a miss for hits beyond the
        // maximum range. This way the free space up to the maximum range will
        // be updated.
        accumulated_range_data_.misses.push_back(
            range_data_in_first_tracking.origin +
            options_.max_range() / range * delta);
      }
    }
  }
  ++num_accumulated_;

  if (num_accumulated_ >= options_.scans_per_accumulation()) {
    num_accumulated_ = 0;
    return AddAccumulatedRangeData(
        time, sensor::TransformRangeData(accumulated_range_data_,
            tracking_delta.inverse()));
  }
  return nullptr;
}

std::unique_ptr<LocalTrajectoryBuilder::InsertionResult> ContinuousLocalTrajectoryBuilder::AddAccumulatedRangeData(
    common::Time time, const sensor::RangeData& range_data_in_tracking) {
  if (imu_tracker_ == nullptr) {
    LOG(INFO)<< "ImuTracker not yet initialized.";
    return nullptr;
  }
  LOG(INFO)<<"acc called";
  const sensor::RangeData filtered_range_data = {
    range_data_in_tracking.origin,
    sensor::VoxelFiltered(range_data_in_tracking.returns,
        options_.voxel_filter_size()),
    sensor::VoxelFiltered(range_data_in_tracking.misses,
        options_.voxel_filter_size())};

  if (first_scan_){
    LOG(INFO)<<"first scan";
    scan_matcher_submap_ =
    active_submaps_.submaps().front();
    first_scan_ = false;
  }
  scan_matcher_submap_ =
  active_submaps_.submaps().front();
  transform::Rigid3d initial_pose = scan_matcher_submap_->local_pose().inverse()*pose_estimate_;
  sensor::AdaptiveVoxelFilter adaptive_voxel_filter(
      options_.high_resolution_adaptive_voxel_filter_options());
  const sensor::PointCloud filtered_point_cloud_in_tracking =
  adaptive_voxel_filter.Filter(filtered_range_data.returns);
  transform::Rigid3d matched_pose = scan_matcher_submap_->local_pose().inverse() * pose_estimate_;
//  real_time_correlative_scan_matcher_->Match(
//          initial_pose, filtered_point_cloud_in_tracking,
//          matching_submap->high_resolution_hybrid_grid(), &matched_pose);

  sensor::AdaptiveVoxelFilter low_resolution_adaptive_voxel_filter(
      options_.low_resolution_adaptive_voxel_filter_options());
  const sensor::PointCloud low_resolution_point_cloud_in_tracking =
  low_resolution_adaptive_voxel_filter.Filter(filtered_range_data.returns);

  transform::Rigid3d last_pose = pose_estimate_;
  transform::Rigid3d pose_observation_in_submap;
  ceres::Solver::Summary summary;
  ceres_scan_matcher_->Match(
      scan_matcher_submap_->local_pose().inverse() * pose_estimate_,
      matched_pose,
      {{&filtered_point_cloud_in_tracking,
        &scan_matcher_submap_->high_resolution_hybrid_grid()},
       {&low_resolution_point_cloud_in_tracking,
        &scan_matcher_submap_->low_resolution_hybrid_grid()}},
      &pose_observation_in_submap, &summary);
//  real_time_correlative_scan_matcher_->Match(
//      scan_matcher_submap_->local_pose().inverse() * pose_estimate_, filtered_point_cloud_in_tracking,
//      scan_matcher_submap_->high_resolution_hybrid_grid(), &pose_observation_in_submap);
  pose_estimate_ = scan_matcher_submap_->local_pose() * pose_observation_in_submap;

  if (last_scan_match_time_ > common::Time::min() &&
      time > last_scan_match_time_) {
    const double delta_t = common::ToSeconds(time - last_scan_match_time_);
    // This adds the observed difference in velocity that would have reduced the
    // error to zero.
//      LOG(INFO)<<"update vel";
    velocity_estimate_ +=
    (pose_estimate_.translation() - last_pose.translation()) /
    delta_t;
  }
  last_scan_match_time_ = time;
//  LOG(INFO)<<"matched pose: "<<matched_pose;
  last_pose_estimate_ = {
    time, pose_estimate_,
    sensor::TransformPointCloud(filtered_range_data.returns,
        pose_estimate_.cast<float>())};

  //insert control points
//  PoseEstimate control_point = {time, pose_estimate_, {}};
//  pose_vec_.push_back(control_point);
  last_control_point_time_ = time;
  if (pose_vec_.size() == 5)
  {
    LOG(INFO)<<"create spline";
    sensor::RangeData complete_scan;
    for (PoseAndRangeData& pose_and_range : local_pose_graph_.createSplineFromControlVector(pose_vec_, range_data_vector_, &scan_matcher_submap_->high_resolution_hybrid_grid(), ceres_scan_matcher_)) {
      sensor::RangeData range_data = sensor::TransformRangeData(pose_and_range.range_data,
          pose_and_range.pose.cast<float>());
      for (const Eigen::Vector3f& hit : range_data.returns) {
        const Eigen::Vector3f delta = hit - range_data.origin;
        const float range = delta.norm();
        if (range >= options_.min_range()) {
          if (range <= options_.max_range()) {
            complete_scan.returns.push_back(hit);
          } else {
            // We insert a ray cropped to 'max_range' as a miss for hits beyond the
            // maximum range. This way the free space up to the maximum range will
            // be updated.
            complete_scan.misses.push_back(
                range_data.origin +
                options_.max_range() / range * delta);
          }
        }
      }
    }
    pose_vec_.clear();
    range_data_vector_.clear();
    // Querying the active submaps must be done here before calling
    // InsertRangeData() since the queried values are valid for next insertion.
    std::vector<std::shared_ptr<const Submap>> insertion_submaps;
    for (std::shared_ptr<Submap> submap : active_submaps_.submaps()) {
      insertion_submaps.push_back(submap);
    }
    active_submaps_.InsertRangeData(
        complete_scan,
        imu_tracker_->orientation());
    scan_matcher_submap_ = active_submaps_.submaps().front();
    return std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
        new LocalTrajectoryBuilder::InsertionResult {time, range_data_in_tracking, pose_estimate_,
          std::move(insertion_submaps)});
  }
//  return InsertIntoSubmap(time, filtered_range_data, pose_estimate_);
  LOG(INFO)<<"finished";
//  return InsertIntoSubmap(time, filtered_range_data, pose_estimate_);
  //return std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(new LocalTrajectoryBuilder::InsertionResult{common::Time(),ranges, PoseEstimate() ,});
  //return std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(new LocalTrajectoryBuilder::InsertionResult());
  return nullptr;

}
std::unique_ptr<LocalTrajectoryBuilder::InsertionResult> ContinuousLocalTrajectoryBuilder::InsertIntoSubmap(
    common::Time time, const sensor::RangeData& range_data_in_tracking,
    const transform::Rigid3d& pose_observation) {
  // Querying the active submaps must be done here before calling
  // InsertRangeData() since the queried values are valid for next insertion.
  std::vector<std::shared_ptr<const Submap>> insertion_submaps;
  for (std::shared_ptr<Submap> submap : active_submaps_.submaps()) {
    insertion_submaps.push_back(submap);
  }
  active_submaps_.InsertRangeData(
      sensor::TransformRangeData(range_data_in_tracking,
                                 pose_observation.cast<float>()),
      imu_tracker_->orientation());
  return std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
      new LocalTrajectoryBuilder::InsertionResult { time,
          range_data_in_tracking, pose_observation, std::move(insertion_submaps) });
}

const mapping::GlobalTrajectoryBuilderInterface::PoseEstimate& ContinuousLocalTrajectoryBuilder::pose_estimate() const {
  return last_pose_estimate_;
}
}
}
