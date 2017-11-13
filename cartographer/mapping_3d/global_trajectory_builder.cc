/*
 * Copyright 2016 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "cartographer/mapping_3d/global_trajectory_builder.h"

namespace cartographer {
namespace mapping_3d {

GlobalTrajectoryBuilder::GlobalTrajectoryBuilder(
    const proto::LocalTrajectoryBuilderOptions& options,
    const int trajectory_id, SparsePoseGraph* sparse_pose_graph)
    : trajectory_id_(trajectory_id),
      sparse_pose_graph_(sparse_pose_graph),
      local_trajectory_builder_(options) {}

GlobalTrajectoryBuilder::~GlobalTrajectoryBuilder() {}

void GlobalTrajectoryBuilder::AddImuData(
    const common::Time time, const Eigen::Vector3d& linear_acceleration,
    const Eigen::Vector3d& angular_velocity) {
  local_trajectory_builder_.AddImuData(time, linear_acceleration,
                                       angular_velocity);
  sparse_pose_graph_->AddImuData(trajectory_id_, time, linear_acceleration,
                                 angular_velocity);
}


/*Selects which local Trajectory Interface should be used*/
void GlobalTrajectoryBuilder::AddRangefinderData(
    const common::Time time, const Eigen::Vector3f& origin,
    const sensor::PointCloud& ranges) {

  /*Continuous Local Trajectory Builder Result*/
//  std::vector<std::unique_ptr<LocalTrajectoryBuilder::InsertionResult>> insertion_result =
//      local_trajectory_builder_.AddRangefinderData(time, origin, ranges);
//  if (insertion_result.size() == 0) {
//    return;
//  }
//  for (auto& result : insertion_result)
//  {
////    LOG(INFO)<<"new result: "<<"time: "<<result->time<<" pose: "<<result->pose_observation
////        <<" range data size: "<<result->range_data_in_tracking.returns.size();
//    sparse_pose_graph_->AddScan(
//        result->time, result->range_data_in_tracking,
//        result->pose_observation, trajectory_id_,
//        std::move(result->insertion_submaps));
//  }
  /*Local Trajectory Builder Result*/
  auto insertion_result =
      local_trajectory_builder_.AddRangefinderData(time, origin, ranges);
  if (insertion_result == nullptr) {
    return;
  }
  sparse_pose_graph_->AddScan(
      insertion_result->time, insertion_result->range_data_in_tracking,
      insertion_result->pose_observation, trajectory_id_,
      std::move(insertion_result->insertion_submaps));
}

void GlobalTrajectoryBuilder::AddOdometerData(const common::Time time,
                                              const transform::Rigid3d& pose) {
  local_trajectory_builder_.AddOdometerData(time, pose);
}

void GlobalTrajectoryBuilder::AddPlaneData(const common::Time time,
					   const Eigen::Vector4d& coefficients)
{
  local_trajectory_builder_.AddPlaneData(time, coefficients);
}

const GlobalTrajectoryBuilder::PoseEstimate&
GlobalTrajectoryBuilder::pose_estimate() const {
  return local_trajectory_builder_.pose_estimate();
}

}  // namespace mapping_3d
}  // namespace cartographer
