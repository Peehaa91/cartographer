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

#include "cartographer/mapping/imu_tracker.h"

#include <cmath>
#include <limits>

#include "cartographer/common/math.h"
#include "cartographer/transform/transform.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

ImuTracker::ImuTracker(const double imu_gravity_time_constant,
                       const common::Time time)
    : imu_gravity_time_constant_(imu_gravity_time_constant),
      time_(time),
      last_linear_acceleration_time_(common::Time::min()),
      orientation_(Eigen::Quaterniond::Identity()),
      gravity_vector_(Eigen::Vector3d::UnitZ()),
      linear_acceleration_correction_(Eigen::Vector3d(0, 0 , 0)),
      first_gravity_(Eigen::Vector3d(0,0, 9.81)),
      imu_angular_velocity_(Eigen::Vector3d::Zero()) {

//    double roll = (M_PI- 3.00584)/2.0;
//    double pitch = (M_PI - 3.12515)/2.0;
//    double yaw = (M_PI - -3.12515)/2.0;
//  orientation_ = Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX())
//  * Eigen::AngleAxisd(pitch,  Eigen::Vector3d::UnitY())
//  * Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ());
//////  orientation_ = orientation_.inverse();
//  gravity_vector_ = orientation_.inverse() * first_gravity_;
}

void ImuTracker::Advance(const common::Time time) {
  CHECK_LE(time_, time);
  const double delta_t = common::ToSeconds(time - time_);
  const Eigen::Quaterniond rotation =
      transform::AngleAxisVectorToRotationQuaternion(
          Eigen::Vector3d(imu_angular_velocity_ * delta_t));
  orientation_ = (orientation_ * rotation).normalized();
  gravity_vector_ = rotation.inverse() * gravity_vector_;
 // linear_acceleration_ = rotation.inverse() * linear_acceleration_;
  time_ = time;
}

void ImuTracker::AddImuLinearAccelerationObservation(
    const Eigen::Vector3d& imu_linear_acceleration) {
  // Update the 'gravity_vector_' with an exponential moving average using the
  // 'imu_gravity_time_constant'.
//  if ( last_linear_acceleration_time_ > common::Time::min())
//  {
  const double delta_t =
      last_linear_acceleration_time_ > common::Time::min()
          ? common::ToSeconds(time_ - last_linear_acceleration_time_)
          : std::numeric_limits<double>::infinity();
//  if (delta_t < std::numeric_limits<double>::infinity())
//    linear_acceleration_ = imu_linear_acceleration - gravity_vector_;
  last_linear_acceleration_time_ = time_;
  const double alpha = 1. - std::exp(-delta_t / imu_gravity_time_constant_);
  gravity_vector_ =
      (1. - alpha) * gravity_vector_ + alpha * imu_linear_acceleration;
//  }
//  else
//    last_linear_acceleration_time_ = time_;



  // Change the 'orientation_' so that it agrees with the current
  // 'gravity_vector_'.

  const Eigen::Quaterniond rotation = Eigen::Quaterniond::FromTwoVectors(
      gravity_vector_, orientation_.inverse() * Eigen::Vector3d::UnitZ());
//  linear_acceleration_correction_ = alpha * imu_linear_acceleration;
  orientation_ = (orientation_ * rotation).normalized();
//  if (delta_t == std::numeric_limits<double>::infinity()){
//    first_gravity_ = Eigen::Vector3d(0, 0, 9.81);
//    LOG(INFO)<<orientation_ * first_gravity_;
//    linear_acceleration_ = imu_linear_acceleration - orientation_.inverse() * first_gravity_;
//    LOG(INFO)<<"first acc: "<<linear_acceleration_;
//  }
//  else{
//    linear_acceleration_ = imu_linear_acceleration - orientation_.inverse() * first_gravity_;
//    LOG(INFO)<<"acc: "<<linear_acceleration_;
//  }
//  linear_acceleration_ = imu_linear_acceleration - orientation_.inverse() * first_gravity_;
  CHECK_GT((orientation_ * gravity_vector_).z(), 0.);
  CHECK_GT((orientation_ * gravity_vector_).normalized().z(), 0.99);
}

void ImuTracker::AddImuAngularVelocityObservation(
    const Eigen::Vector3d& imu_angular_velocity) {
  imu_angular_velocity_ = imu_angular_velocity;
}

}  // namespace mapping
}  // namespace cartographer
