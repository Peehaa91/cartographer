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
  for (int i  = 0; i < 5; i++){
    last_pose_vec_.push_back(last_pose_estimate_);
  }
  first_spline_ = true;
}

double counter = 0;
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
  if (range_data_vector_.size() > 0) {
    imu_angular_vel_data_vector_.push_back(
        std::make_pair(time, angular_velocity));
    linear_acc_data_vector_.push_back(
        std::make_pair(time, linear_acceleration));
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
    for (PoseEstimate& control_point : last_pose_vec_){
      Eigen::Quaterniond control_point_orientation = control_point.pose.rotation();
      control_point_orientation = control_point_orientation*
          last_orientation.inverse() * imu_tracker_->orientation();
      control_point.pose =  transform::Rigid3d(control_point.pose.translation(),
                                               control_point_orientation);
    }
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

std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>> ContinuousLocalTrajectoryBuilder::AddRangefinderData(
    const common::Time time, const Eigen::Vector3f& origin,
    const sensor::PointCloud& ranges) {
  if (imu_tracker_ == nullptr) {
    LOG(INFO)<< "ImuTracker not yet initialized.";
    return std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>>();
  }
  Predict(time);
  if (first_scan_) {
    LOG(INFO)<<"first scan";
//    scan_matcher_submap_ =
//    active_submaps_.submaps().front();
//    PoseEstimate control_point = {time, scan_matcher_submap_->local_pose().inverse()*pose_estimate_, {}};
    PoseEstimate control_point = {time, pose_estimate_, {}};
    pose_vec_.push_back(control_point);
    first_scan_ = false;
  }
//  if ((range_data_vector_.empty() && pose_vec_.empty())){// || control_point_counter_ == 32){
//    pose_vec_.push_back( {time, pose_estimate_, {}});
//    control_point_counter_ = 0;
//  }
//  LOG(INFO)<<"before:"<<ranges.size();
  const sensor::PointCloud filtered_point_cloud = sensor::VoxelFiltered(ranges,
      options_.voxel_filter_size());
  sensor::AdaptiveVoxelFilter adaptive_voxel_filter(
      options_.high_resolution_adaptive_voxel_filter_options());
  sensor::PointCloud filtered_point_cloud_in_tracking_first = sensor::VoxelFiltered(filtered_point_cloud,
                        options_.voxel_filter_size());
  const sensor::PointCloud filtered_point_cloud_in_tracking =
  adaptive_voxel_filter.Filter(filtered_point_cloud_in_tracking_first);
//  LOG(INFO)<<"after:"<<filtered_point_cloud_in_tracking.size();
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
  if (num_accumulated_ >= 320 || (first_spline_ && (num_accumulated_ >= 320))) {
    num_accumulated_ = 0;
    return AddAccumulatedRangeData(
        time, sensor::TransformRangeData(accumulated_range_data_,
            tracking_delta.inverse()));
  }
  return std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>>();
}

std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>> ContinuousLocalTrajectoryBuilder::AddAccumulatedRangeData(
    common::Time time, const sensor::RangeData& range_data_in_tracking) {
  if (imu_tracker_ == nullptr) {
    LOG(INFO)<< "ImuTracker not yet initialized.";
    return std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>>();
  }
  const sensor::RangeData filtered_range_data = {
    range_data_in_tracking.origin,
    sensor::VoxelFiltered(range_data_in_tracking.returns,
        options_.voxel_filter_size()),
    sensor::VoxelFiltered(range_data_in_tracking.misses,
        options_.voxel_filter_size())};

  std::shared_ptr<const SubmapDecay> matching_submap =
  active_submaps_.submaps().front();
  LOG(INFO)<<"submap size:"<<active_submaps_.submaps().size();
  transform::Rigid3d initial_pose = matching_submap->local_pose().inverse()*pose_estimate_;
  transform::Rigid3d last_pose = pose_estimate_;

  //insert control points
  PoseEstimate control_point;
  if (!first_spline_)
  {
    //Propogate Control Points
    if (pose_vec_.size() == 1)
      pose_vec_[0].pose = matching_submap->local_pose().inverse() * pose_vec_[0].pose;
    Eigen::Vector3d translation = last_pose_vec_[pose_vec_.size()].pose.translation()
        + control_point_vel_[pose_vec_.size()] * common::ToSeconds(time - last_pose_vec_[pose_vec_.size()].time);
//    LOG(INFO)<<pose_vec_.size();
    transform::Rigid3d new_pose = transform::Rigid3d(translation, pose_estimate_.rotation());
//    transform::Rigid3d new_pose = transform::Rigid3d(last_pose_vec_[last_pose_vec_.size() -1].pose.translation(),
//                                                     pose_estimate_.rotation());
//    transform::Rigid3d new_pose = pose_estimate_;
//    transform::Rigid3d new_pose = transform::Rigid3d(pose_estimate_.translation(),last_pose_vec_[pose_vec_.size()].pose.rotation());
    control_point = {time, matching_submap->local_pose().inverse()*new_pose, {}};
//    LOG(INFO)<<"delta trans:"<<control_point_vel_[pose_vec_.size()] * common::ToSeconds(time - last_pose_vec_[pose_vec_.size()].time);
//    LOG(INFO)<<"vel:"<<control_point_vel_[pose_vec_.size()];
//    LOG(INFO)<<"delta time:"<<common::ToSeconds(time - last_pose_vec_[pose_vec_.size()].time);
//    LOG(INFO)<<"new pose:"<<new_pose.translation();
//    LOG(INFO)<<"new control:"<<control_point.pose;
//    LOG(INFO)<<"local pose:"<<scan_matcher_submap_->local_pose();
//    LOG(INFO)<<"id_index:"<<active_submaps_.matching_index();
    pose_vec_.push_back(control_point);
  }
  else{
    control_point = {time, matching_submap->local_pose().inverse()*pose_estimate_, {}};
    pose_vec_.push_back(control_point);
  }

  last_control_point_time_ = time;
  const HybridDecayGrid* decay_grid_high = &matching_submap->high_resolution_hybrid_decay_grid();
  const HybridDecayGrid* decay_grid_low = &matching_submap->low_resolution_hybrid_decay_grid();
  std::vector<const HybridDecayGrid*> decay_grids = {decay_grid_high, decay_grid_low};
  std::vector<const HybridGrid*> hybrid_grids = {&matching_submap->high_resolution_hybrid_grid(),
    &matching_submap->low_resolution_hybrid_grid()};
  sensor::AdaptiveVoxelFilter adaptive_voxel_filter(
      options_.high_resolution_adaptive_voxel_filter_options());
  const sensor::PointCloud filtered_point_cloud_in_tracking =
      adaptive_voxel_filter.Filter(range_data_in_tracking.returns);


  if (first_spline_){
    LOG(INFO)<<"first spline";
    sensor::RangeData complete_scan = filtered_range_data;
    first_spline_ = false;
    last_pose_vec_ = pose_vec_;
    pose_vec_.clear();
    range_data_vector_.clear();
    imu_angular_vel_data_vector_.clear();
    linear_acc_data_vector_.clear();
    last_scan_match_time_ = time;
    pose_vec_.push_back(control_point);
//    last_pose_estimate_ = {
//      time, pose_estimate_,
//      sensor::TransformPointCloud(complete_scan.returns,
//          pose_estimate_.cast<float>())};
    last_pose_estimate_ = {
        time, pose_estimate_,
        complete_scan.returns};
    std::vector<std::shared_ptr<const Submap>> insertion_submaps;
    for (std::shared_ptr<Submap> submap : active_submaps_.submaps()) {
      insertion_submaps.push_back(submap);
    }
    active_submaps_.InsertRangeData(sensor::TransformRangeData(
        complete_scan, pose_estimate_.cast<float>()),
        pose_estimate_.rotation(),
        true);
    std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>> result;
    result.push_back(std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
        new LocalTrajectoryBuilder::InsertionResult {time, complete_scan, pose_estimate_,
          std::move(insertion_submaps)}));
    return result;
  }
  if (pose_vec_.size() == 5)
  {
    LOG(INFO)<<"create spline";
    sensor::RangeData complete_scan;
    complete_scan.origin = pose_estimate_.translation().cast<float>();
    for (auto& acc : linear_acc_data_vector_)
    {
      acc = std::make_pair(acc.first, matching_submap->local_pose().rotation().inverse()* acc.second);
    }
    std::vector<PoseAndRangeData>pose_and_range_vec = local_pose_graph_.createSplineFromControlVector(pose_vec_,
        range_data_vector_,
        imu_angular_vel_data_vector_,
        linear_acc_data_vector_,
        hybrid_grids,
        ceres_scan_matcher_,
        decay_grids);
    for (int i = 0; i < 5; i++)
    {
      pose_vec_[i].pose = matching_submap->local_pose() * pose_vec_[i].pose;
      control_point_vel_[i] = (pose_vec_[i].pose.translation() - last_pose_vec_[i].pose.translation())/
          common::ToSeconds(pose_vec_[i].time - last_pose_vec_[i].time);
      LOG(INFO)<<"vel "<<i<<": "<<control_point_vel_[i];
      last_pose_vec_[i].pose = matching_submap->local_pose()*last_pose_vec_[i].pose;
    }
    last_pose_vec_ = pose_vec_;


    //Complete Scan for last_pose_estimate Interface
    transform::Rigid3d last_pose;
    common::Time last_time = common::Time::min();
    for (PoseAndRangeData& pose_and_range : pose_and_range_vec) {
      //Transform range data to map
      sensor::RangeData range_data = sensor::TransformRangeData(pose_and_range.range_data,
          matching_submap->local_pose().cast<float>()*pose_and_range.pose.cast<float>());
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
//    velocity_estimate_ = (pose_vec_[4].pose.translation() - pose_vec_[3].pose.translation())/
//        common::ToSeconds(pose_vec_[4].time - pose_vec_[3].time);
//    LOG(INFO)<<"new vel: "<<velocity_estimate_;
    pose_estimate_ = matching_submap->local_pose() * pose_and_range_vec.back().pose;
    pose_vec_.clear();
    range_data_vector_.clear();
    imu_angular_vel_data_vector_.clear();
    linear_acc_data_vector_.clear();

    last_scan_match_time_ = time;
    LOG(INFO)<<"old pose:"<<pose_estimate_;

    LOG(INFO)<<"new pose:"<<pose_estimate_;
    control_point = {time, pose_estimate_, {}};
    pose_vec_.push_back(control_point);
//    last_pose_estimate_ = {
//      time, pose_estimate_,
//      sensor::TransformPointCloud(complete_scan.returns,
//          pose_estimate_.cast<float>())};
    last_pose_estimate_ = {
        time, pose_estimate_,
        complete_scan.returns};
//
//    active_submaps_.InsertRangeData(
//        complete_scan,
//        pose_estimate_.rotation());


    // mapping and connection for the global slam result for the single "scans"
    bool last_data = false;
    int range_data_counter = 0;
    std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>> insertion_result_vec;
    for (PoseAndRangeData& pose_and_range : pose_and_range_vec) {
      //Transform range data to map
//      sensor::RangeData range_data = sensor::TransformRangeData(pose_and_range.range_data,
//                                                                matching_submap->local_pose().cast<float>()*
//                                                                pose_and_range.pose.cast<float>());
      sensor::RangeData range_data = pose_and_range.range_data;
      sensor::RangeData temp_cloud;
      temp_cloud.origin = range_data.origin;
      for (const Eigen::Vector3f& hit : range_data.returns) {
        const Eigen::Vector3f delta = hit - range_data.origin;
        const float range = delta.norm();
        if (range >= options_.min_range()) {
          if (range <= options_.max_range()) {
            temp_cloud.returns.push_back(hit);
          } else {
            // We insert a ray cropped to 'max_range' as a miss for hits beyond the
            // maximum range. This way the free space up to the maximum range will
            // be updated.
            temp_cloud.misses.push_back(
                range_data.origin +
                options_.max_range() / range * delta);
          }
        }
      }
      if (range_data_counter == 0)
        last_data = true;

      range_data_counter++;
      transform::Rigid3d transform = matching_submap->local_pose()*pose_and_range.pose;
      active_submaps_.InsertRangeData(sensor::TransformRangeData(temp_cloud,transform.cast<float>()),
          imu_tracker_->orientation(),
          last_data);
      last_data = false;
      std::vector<std::shared_ptr<const Submap>> insertion_submaps;
      for (std::shared_ptr<Submap> submap : active_submaps_.submaps()) {
        insertion_submaps.push_back(submap);
      }
      temp_cloud.origin = range_data.origin;
      insertion_result_vec.push_back(std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
              new LocalTrajectoryBuilder::InsertionResult {pose_and_range.time, temp_cloud, matching_submap->local_pose()*pose_and_range.pose,
                std::move(insertion_submaps)}));
//      LOG(INFO)<<matching_submap->local_pose();
    }
//    std::vector<std::shared_ptr<const Submap>> insertion_submaps;
//    for (std::shared_ptr<Submap> submap : active_submaps_.submaps()) {
//      insertion_submaps.push_back(submap);
//    }
//    insertion_result_vec.push_back(std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
//        new LocalTrajectoryBuilder::InsertionResult {time, complete_scan, pose_estimate_,
//          std::move(insertion_submaps)}));
    return insertion_result_vec;
    //end insertion vec
//    return std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
//        new LocalTrajectoryBuilder::InsertionResult {time, complete_scan, pose_estimate_,
//          std::move(insertion_submaps)});
  }
  return std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>>();

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
//  active_submaps_.InsertRangeData(
//      sensor::TransformRangeData(range_data_in_tracking,
//                                 pose_observation.cast<float>()),
//      imu_tracker_->orientation());
  return std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>(
      new LocalTrajectoryBuilder::InsertionResult { time,
          range_data_in_tracking, pose_observation, std::move(insertion_submaps) });
}

const mapping::GlobalTrajectoryBuilderInterface::PoseEstimate& ContinuousLocalTrajectoryBuilder::pose_estimate() const {
  return last_pose_estimate_;
}
}
}
