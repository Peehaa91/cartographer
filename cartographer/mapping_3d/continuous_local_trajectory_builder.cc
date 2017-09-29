#include <cartographer/mapping_3d/continuous_local_trajectory_builder.h>
#include "cartographer/mapping_3d/scan_matching/proto/ceres_scan_matcher_options.pb.h"


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
  if (range_data_vector_.size() > 0){
    imu_angular_vel_data_vector_.push_back(std::make_pair(time, angular_velocity));
    linear_acc_data_vector_.push_back(std::make_pair(time, linear_acceleration));
  }
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
  CHECK_LE(last_update_time_, time);
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
    acceleration_estimate_ = imu_tracker_->linear_acceleration();
//    velocity_estimate_ += acceleration_estimate_ * delta_t;
//    velocity_estimate_ = rotation.inverse() * velocity_estimate_;
//      last_pose_estimate_ = PoseEstimate(time, pose,  {});
    pose_estimate_ = pose;
//    vel_est_ += imu_tracker_->linear_acceleration() * delta_t;
//    LOG(INFO)<<"vel_est: "<<vel_est_;
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
  if (first_scan_){
     LOG(INFO)<<"first scan";
     scan_matcher_submap_ =
     active_submaps_.submaps().front();
     PoseEstimate control_point = {time, scan_matcher_submap_->local_pose().inverse()*pose_estimate_, {}};
     pose_vec_.push_back(control_point);
     first_scan_ = false;
   }
//  if ((range_data_vector_.empty() && pose_vec_.empty())){// || control_point_counter_ == 32){
//    pose_vec_.push_back( {time, pose_estimate_, {}});
//    control_point_counter_ = 0;
//  }
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
  const sensor::RangeData filtered_range_data = {
    range_data_in_tracking.origin,
    sensor::VoxelFiltered(range_data_in_tracking.returns,
        options_.voxel_filter_size()),
    sensor::VoxelFiltered(range_data_in_tracking.misses,
        options_.voxel_filter_size())};

  scan_matcher_submap_ =
  active_submaps_.submaps().front();
  transform::Rigid3d initial_pose = scan_matcher_submap_->local_pose().inverse()*pose_estimate_;
  transform::Rigid3d last_pose = pose_estimate_;

  //insert control points
  PoseEstimate control_point = {time, scan_matcher_submap_->local_pose().inverse()*pose_estimate_, {}};
  pose_vec_.push_back(control_point);
  last_control_point_time_ = time;
  const HybridDecayGrid* decay_grid = &scan_matcher_submap_->high_resolution_hybrid_decay_grid();

  if (pose_vec_.size() == 5)
  {
    LOG(INFO)<<"create spline";
    sensor::RangeData complete_scan;
    std::vector<PoseAndRangeData>pose_and_range_vec = local_pose_graph_.createSplineFromControlVector(pose_vec_,
                                                                                                      range_data_vector_,
                                                                                                      imu_angular_vel_data_vector_,
                                                                                                      linear_acc_data_vector_,
                                                                                                      &scan_matcher_submap_->high_resolution_hybrid_grid(),
                                                                                                      ceres_scan_matcher_,
                                                                                                      decay_grid);
    transform::Rigid3d last_pose;
    common::Time last_time = common::Time::min();
    for (PoseAndRangeData& pose_and_range : pose_and_range_vec) {
      sensor::RangeData range_data = sensor::TransformRangeData(pose_and_range.range_data,
                                                                scan_matcher_submap_->local_pose().cast<float>()*pose_and_range.pose.cast<float>());
//      if (last_time > common::Time::min()){
//        double delta_t = common::ToSeconds(pose_and_range.time - last_time);
//        velocity_estimate_ +=
//        (pose_and_range.pose.translation()
//            - (last_pose).translation()) /
//            (delta_t);
//        LOG(INFO)<<"new vel: "<<velocity_estimate_<<"pose: "<<pose_and_range.pose.translation()<<" last pose: "<<(last_pose).translation()<<"delta_t "<<delta_t;
//      }
      last_time = pose_and_range.time;
      last_pose = pose_and_range.pose;
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
    imu_angular_vel_data_vector_.clear();
    linear_acc_data_vector_.clear();
    // Querying the active submaps must be done here before calling
    // InsertRangeData() since the queried values are valid for next insertion.
    if (last_scan_match_time_ > common::Time::min() && time > last_scan_match_time_) {
      const double delta_t = common::ToSeconds(time - pose_and_range_vec[pose_and_range_vec.size() - 2].time);
//      const double delta_t = common::ToSeconds(pose_and_range_vec.back().time
//              - pose_and_range_vec[pose_and_range_vec.size()-2].time);
//      velocity_estimate_ =
//      ((scan_matcher_submap_->local_pose()*pose_and_range_vec.back().pose).translation()
//          - scan_matcher_submap_->local_pose()* pose_and_range_vec[pose_and_range_vec.size()-2].pose.translation()) /
//          (delta_t);
//      LOG(INFO)<<"vel before: "<<velocity_estimate_;
//      velocity_estimate_ =
//      ((scan_matcher_submap_->local_pose()*pose_and_range_vec.back().pose).translation()
//          - (scan_matcher_submap_->local_pose()* pose_and_range_vec[pose_and_range_vec.size()-2].pose).translation()) /
//          (delta_t);
      LOG(INFO)<<"new vel: "<<velocity_estimate_;
//      velocity_estimate_ = Eigen::Vector3d::Zero();
    }
    last_scan_match_time_ = time;
    LOG(INFO)<<"old pose:"<<pose_estimate_;
    pose_estimate_ = scan_matcher_submap_->local_pose() * pose_and_range_vec.back().pose;
    LOG(INFO)<<"new pose:"<<pose_estimate_;
    control_point = {time, scan_matcher_submap_->local_pose().inverse()*pose_estimate_, {}};
    pose_vec_.push_back(control_point);
    last_pose_estimate_ = {
      time, pose_estimate_,
      sensor::TransformPointCloud(filtered_range_data.returns,
          pose_estimate_.cast<float>())};
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
//  LOG(INFO)<<"finished";
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
