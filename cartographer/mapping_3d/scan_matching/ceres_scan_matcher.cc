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

#include "cartographer/mapping_3d/scan_matching/ceres_scan_matcher.h"

#include <string>
#include <utility>
#include <vector>

#include "cartographer/common/ceres_solver_options.h"
#include "cartographer/common/make_unique.h"
#include "cartographer/mapping_3d/ceres_pose.h"
#include "cartographer/mapping_3d/scan_matching/occupied_space_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/rotation_delta_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/translation_delta_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/knot_translation_delta_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_weight_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_control_points_occ_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_control_point_auto_cost_functor.h"
//#include "cartographer/mapping_3d/scan_matching/nurbs_control_point_dynamic_auto_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_knot_weights_auto_cost_functor.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_rotation_vel_delta_functor.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_translation_acceleration_delta_functor.h"
#include "cartographer/mapping_3d/scan_matching/free_space_cost_functor.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/transform/transform.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "cartographer/mapping_3d/ray_tracer.h"
#include "cartographer/mapping_3d/scan_matching/nurbs_control_points_free_space_cost_functor.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {
namespace {

struct YawOnlyQuaternionPlus {
  template <typename T>
  bool operator()(const T* x, const T* delta, T* x_plus_delta) const {
    const T clamped_delta = common::Clamp(delta[0], T(-0.5), T(0.5));
    T q_delta[4];
    q_delta[0] = ceres::sqrt(1. - clamped_delta * clamped_delta);
    q_delta[1] = T(0.);
    q_delta[2] = T(0.);
    q_delta[3] = clamped_delta;
    ceres::QuaternionProduct(q_delta, x, x_plus_delta);
    return true;
  }
};

}  // namespace

proto::CeresScanMatcherOptions CreateCeresScanMatcherOptions(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::CeresScanMatcherOptions options;
  for (int i = 0;; ++i) {
    const string lua_identifier = "occupied_space_weight_" + std::to_string(i);
    if (!parameter_dictionary->HasKey(lua_identifier)) {
      break;
    }
    options.add_occupied_space_weight(
        parameter_dictionary->GetDouble(lua_identifier));
  }
  options.set_translation_weight(
      parameter_dictionary->GetDouble("translation_weight"));
  options.set_rotation_weight(
      parameter_dictionary->GetDouble("rotation_weight"));
  options.set_only_optimize_yaw(
      parameter_dictionary->GetBool("only_optimize_yaw"));
  LOG(INFO)<<"has ray tracer: "<<options.has_ray_tracer_line_size();
  options.set_ray_tracer_line_size(
      parameter_dictionary->GetInt("ray_tracer_line_size"));
  *options.mutable_ceres_solver_options() =
      common::CreateCeresSolverOptionsProto(
          parameter_dictionary->GetDictionary("ceres_solver_options").get());
  return options;
}

CeresScanMatcher::CeresScanMatcher(
    const proto::CeresScanMatcherOptions& options)
    : options_(options),
      ceres_solver_options_(
          common::CreateCeresSolverOptions(options.ceres_solver_options())) {
  ceres_solver_options_.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  ceres_solver_options_.num_linear_solver_threads = options_.ceres_solver_options().num_threads();
  ray_tracer_ = std::make_shared<RayTracer>(RayTracer(options.ray_tracer_line_size()));
}
int counter = 0;
double sum = 0.0;
void CeresScanMatcher::Match(const transform::Rigid3d& previous_pose,
                             const transform::Rigid3d& initial_pose_estimate,
                             const std::vector<PointCloudAndHybridGridPointers>&
                                 point_clouds_and_hybrid_grids,
                             transform::Rigid3d* const pose_estimate,
                             ceres::Solver::Summary* const summary) {
  ceres::Problem problem;
  CeresPose ceres_pose(
      initial_pose_estimate, nullptr /* translation_parameterization */,
      options_.only_optimize_yaw()
          ? std::unique_ptr<ceres::LocalParameterization>(
                common::make_unique<ceres::AutoDiffLocalParameterization<
                    YawOnlyQuaternionPlus, 4, 1>>())
          : std::unique_ptr<ceres::LocalParameterization>(
                common::make_unique<ceres::QuaternionParameterization>()),
      &problem);
  CHECK_EQ(options_.occupied_space_weight_size(),
           point_clouds_and_hybrid_grids.size());
  for (size_t i = 0; i != point_clouds_and_hybrid_grids.size() ; ++i) {
    CHECK_GT(options_.occupied_space_weight(i), 0.);
    const sensor::PointCloud& point_cloud =
        *point_clouds_and_hybrid_grids[i].first;
    const HybridGrid& hybrid_grid = *point_clouds_and_hybrid_grids[i].second;
    LOG(INFO)<<"ceres before: "<<ceres_pose.ToRigid();
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<OccupiedSpaceCostFunctor,
                                        ceres::DYNAMIC, 3, 4>(
            new OccupiedSpaceCostFunctor(
                options_.occupied_space_weight(i) /
                    std::sqrt(static_cast<double>(point_cloud.size())),
                point_cloud, hybrid_grid),
            point_cloud.size()),
        nullptr, ceres_pose.translation(), ceres_pose.rotation());
  }
  CHECK_GT(options_.translation_weight(), 0.);
  problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<TranslationDeltaCostFunctor, 3, 3>(
          new TranslationDeltaCostFunctor(options_.translation_weight(),
                                          previous_pose)),
      nullptr, ceres_pose.translation());
  CHECK_GT(options_.rotation_weight(), 0.);
  problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<RotationDeltaCostFunctor, 3, 4>(
          new RotationDeltaCostFunctor(options_.rotation_weight(),
                                       initial_pose_estimate.rotation())),
      nullptr, ceres_pose.rotation());

  ceres::Solve(ceres_solver_options_, &problem, summary);
  LOG(INFO)<<"ceres after: "<<ceres_pose.ToRigid();
  sum += summary->total_time_in_seconds;
  counter++;
  LOG(INFO)<<"mean secs: "<<sum/counter;
  LOG(INFO)<<summary->FullReport();
  *pose_estimate = ceres_pose.ToRigid();
}

void CeresScanMatcher::MatchWithFreeSpace(const transform::Rigid3d& previous_pose,
                             const transform::Rigid3d& initial_pose_estimate,
                             const std::vector<PointCloudAndHybridGridPointers>&
                                 point_clouds_and_hybrid_grids,
                                 std::vector<const HybridDecayGrid*> decay_grid,
                             transform::Rigid3d* const pose_estimate,
                             ceres::Solver::Summary* const summary) {
  ceres::Problem problem;
  CeresPose ceres_pose(
      initial_pose_estimate, nullptr /* translation_parameterization */,
      options_.only_optimize_yaw()
          ? std::unique_ptr<ceres::LocalParameterization>(
                common::make_unique<ceres::AutoDiffLocalParameterization<
                    YawOnlyQuaternionPlus, 4, 1>>())
          : std::unique_ptr<ceres::LocalParameterization>(
                common::make_unique<ceres::QuaternionParameterization>()),
      &problem);
  CHECK_EQ(options_.occupied_space_weight_size(),
           point_clouds_and_hybrid_grids.size());
  for (size_t i = 0; i != point_clouds_and_hybrid_grids.size(); ++i) {
    CHECK_GT(options_.occupied_space_weight(i), 0.);
    const sensor::PointCloud& point_cloud =
        *point_clouds_and_hybrid_grids[i].first;
    const HybridGrid& hybrid_grid = *point_clouds_and_hybrid_grids[i].second;
    LOG(INFO)<<"ceres before: "<<ceres_pose.ToRigid();
    LOG(INFO)<<"Free Space evaluation";
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<FreeSpaceCostFunctor,
                                        ceres::DYNAMIC, 3, 4>(
            new FreeSpaceCostFunctor(
                options_.occupied_space_weight(i) /
                    std::sqrt(static_cast<double>(point_cloud.size())),
                point_cloud, hybrid_grid, decay_grid[i], options_.ray_tracer_line_size()),
            point_cloud.size()),
        nullptr, ceres_pose.translation(), ceres_pose.rotation());
  }
//  const sensor::PointCloud& point_cloud =
//      *point_clouds_and_hybrid_grids[0].first;
//  const HybridGrid& hybrid_grid = *point_clouds_and_hybrid_grids[0].second;
//  problem.AddResidualBlock(
//      new ceres::AutoDiffCostFunction<FreeSpaceCostFunctor,
//                                      ceres::DYNAMIC, 3, 4>(
//          new FreeSpaceCostFunctor(
//              options_.occupied_space_weight(0) /
//                  std::sqrt(static_cast<double>(point_cloud.size())),
//              point_cloud, hybrid_grid, decay_grid, options_.ray_tracer_line_size()),
//          point_cloud.size()),
//      nullptr, ceres_pose.translation(), ceres_pose.rotation());
  CHECK_GT(options_.translation_weight(), 0.);
  problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<TranslationDeltaCostFunctor, 3, 3>(
          new TranslationDeltaCostFunctor(options_.translation_weight(),
                                          previous_pose)),
      nullptr, ceres_pose.translation());
  CHECK_GT(options_.rotation_weight(), 0.);
  problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<RotationDeltaCostFunctor, 3, 4>(
          new RotationDeltaCostFunctor(options_.rotation_weight(),
                                       initial_pose_estimate.rotation())),
      nullptr, ceres_pose.rotation());

  ceres::Solve(ceres_solver_options_, &problem, summary);
  LOG(INFO)<<"ceres after: "<<ceres_pose.ToRigid();
  sum += summary->total_time_in_seconds;
  counter++;
  LOG(INFO)<<"mean secs: "<<sum/counter;
  LOG(INFO)<<summary->FullReport();
  *pose_estimate = ceres_pose.ToRigid();
}
void CeresScanMatcher::Match(const HybridGrid* hybrid_grid,
                             const common::Time& begin,
                             const common::Time& end,
                             const std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
                             Nurbs<double, 1, 7, KnotType::UNIFORM, WeightType::RATIONAL>& nurbs,
                             Nurbs<double, 1, 7, KnotType::UNIFORM, WeightType::RATIONAL>* nurbs_estimate,
                             ceres::Solver::Summary* summary)
{
  ceres::Problem problem;
  int number_of_residuals = 0;
  for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec)
    number_of_residuals += data.second.returns.size();
  NurbsWeights<double, 1, 7, WeightType::RATIONAL>& weights = nurbs.getWeights();
  boost::multi_array<double, 1> weights_vec = weights.getWeights();
  std::array<double, 5> initial_weights;
  for (int i = 0; i < weights_vec.size(); i++) {
    initial_weights[i] = weights_vec[i];
  }
  const HybridGrid& hybrid_grid_copy = *hybrid_grid;
//  problem.AddResidualBlock(
//          new ceres:: NumericDiffCostFunction<NurbsWeightCostFunctor,
//          ceres::FORWARD,ceres::DYNAMIC, 5>(
//              new NurbsWeightCostFunctor(
//                  options_.occupied_space_weight(0)/std::sqrt(static_cast<double>(number_of_residuals)),
//                  hybrid_grid_copy,
//                  begin,
//                  end,
//                  range_data_vec,
//                  &nurbs), ceres::TAKE_OWNERSHIP,
//                  number_of_residuals),
//          nullptr, initial_weights.data());
  LOG(INFO)<<"add residual";
  boost::multi_array<std::array<double, 7>, 1>& points = nurbs.getPoints();
  problem.AddResidualBlock(
          new ceres:: NumericDiffCostFunction<NurbsControlPointOccSpaceCostFunctor,
          ceres::FORWARD,ceres::DYNAMIC, 7, 7, 7, 7, 7>(
              new NurbsControlPointOccSpaceCostFunctor(
                  options_.occupied_space_weight(0)/std::sqrt(static_cast<double>(number_of_residuals)),
                  hybrid_grid_copy,
                  begin,
                  end,
                  range_data_vec,
                  &nurbs), ceres::TAKE_OWNERSHIP,
                  number_of_residuals),
          nullptr, points[0].data(), points[1].data(), points[2].data(), points[3].data(), points[4].data());
//  for (int i = 0; i < 5; i++)
//  {
//    problem.AddResidualBlock(
//        new ceres::AutoDiffCostFunction<KnotTranslationDeltaCostFunctor, 3, 6>(
//            new KnotTranslationDeltaCostFunctor(options_.translation_weight(),
//                                            points[i])),
//        nullptr, points[i].data());
//  }
//  std::vector<CeresPose> knot_vec;
//  for (int i = 0; i < 5; i++){
//    transform::Rigid3d initial_pose = transform::Rigid3d(Eigen::Vector3d(points[i][0], points[i][1], points[i][2]),
//                                                         Eigen::Quaterniond(points[i][6], points[i][3], points[i][4], points[i][5]));
////    CeresPose ceres_pose = CeresPose(
////        initial_pose, nullptr /* translation_parameterization */,
////        options_.only_optimize_yaw()
////            ? std::unique_ptr<ceres::LocalParameterization>(
////                  common::make_unique<ceres::AutoDiffLocalParameterization<
////                      YawOnlyQuaternionPlus, 4, 1>>())
////            : std::unique_ptr<ceres::LocalParameterization>(
////                  common::make_unique<ceres::QuaternionParameterization>()),
////        &problem);
//    knot_vec.push_back(CeresPose(
//        initial_pose, nullptr /* translation_parameterization */,
//        options_.only_optimize_yaw()
//            ? std::unique_ptr<ceres::LocalParameterization>(
//                  common::make_unique<ceres::AutoDiffLocalParameterization<
//                      YawOnlyQuaternionPlus, 4, 1>>())
//            : std::unique_ptr<ceres::LocalParameterization>(
//                  common::make_unique<ceres::QuaternionParameterization>()),
//        &problem));
//  }
//  problem.AddResidualBlock(
//          new ceres:: AutoDiffCostFunction<NurbsKnotAutoCostFunctor,
//          ceres::DYNAMIC, 3,4,3,4,3,4,3,4>(
//              new NurbsKnotAutoCostFunctor(
//                  options_.occupied_space_weight(0)/std::sqrt(static_cast<double>(number_of_residuals)),
//                  hybrid_grid_copy,
//                  begin,
//                  end,
//                  range_data_vec,
//                  &nurbs),
//                  number_of_residuals),
//          nullptr, knot_vec[0].translation(), knot_vec[0].rotation()
//          , knot_vec[1].translation(), knot_vec[1].rotation()
//          , knot_vec[2].translation(), knot_vec[2].rotation()
//          , knot_vec[3].translation(), knot_vec[3].rotation());

  std::vector<int> const_para = {4, 5};
  for (int i = 0; i < 5; i++) {
    problem.SetParameterLowerBound(points[i].data(), 3, -1);
    problem.SetParameterUpperBound(points[i].data(), 3, 1);
//    std::unique_ptr<ceres::LocalParameterization> subset_para =  common::make_unique<ceres::SubsetParameterization>(7, const_para);
//    problem.SetParameterization(points[i].data(),subset_para.release());
    for (int j = 4; j < 7; j++)
    {
      problem.SetParameterLowerBound(points[i].data(), j, -1);
      problem.SetParameterUpperBound(points[i].data(), j, 1);
    }
  }

//  }
  std::unique_ptr<ceres::LocalParameterization> parametrization = common::make_unique<ceres:: HomogeneousVectorParameterization>(5);
//  problem.SetParameterization(initial_weights.data(),parametrization.release());
  LOG(INFO)<<"number of threads:"<<ceres_solver_options_.num_threads;
  LOG(INFO)<<"before";
   for (int i = 0; i < points.size(); i++) {
     Eigen::Quaterniond quat(points[i][6], points[i][3], points[i][4], points[i][5]);
     Eigen::Matrix3d rot_mat = Eigen::Matrix3d(quat);
     Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
     LOG(INFO)<<"x: "<<points[i][0]<<" y: "<<points[i][1]<<"z: "<<points[i][2]<<" roll: "<<angles(0)<<" pitch: "<<angles(1)<<" yaw: "<<angles(2);
     }
  ceres::Solve(ceres_solver_options_, &problem, summary);
  LOG(INFO)<<summary->FullReport();
  boost::multi_array<double, 1>& weights_estimate = nurbs_estimate->getWeights().getWeights();
  LOG(INFO)<<"after";
  for (int i = 0; i < points.size(); i++) {
    Eigen::Quaterniond quat(points[i][6], points[i][3], points[i][4], points[i][5]);
    Eigen::Matrix3d rot_mat = Eigen::Matrix3d(quat);
    Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
    LOG(INFO)<<"x: "<<points[i][0]<<" y: "<<points[i][1]<<"z: "<<points[i][2]<<" roll: "<<angles(0)<<" pitch: "<<angles(1)<<" yaw: "<<angles(2);
    }

}

void CeresScanMatcher::MatchSplineWithFreeSpace(const std::vector<const HybridGrid*>& hybrid_grid,
                                                const std::vector<const HybridDecayGrid*>& decay_grid,
                                                const common::Time& begin,
                                                const common::Time& end,
                                                const std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
                                                std::vector<std::pair<common::Time, Eigen::Vector3d>>& imu_angular_vel_data_vector,
                                                std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector,
                                                std::vector<transform::Rigid3d>& control_points,
                                                std::vector<double>& weight_vec,
                                                ceres::Solver::Summary* summary)
{
  ceres::Problem problem;
  int number_of_residuals = 0;
  for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec)
    number_of_residuals += data.second.returns.size();
  std::vector<CeresPose> ceres_pose_control_points;
  for (int i = 0; i < control_points.size(); i++){
    if (options_.only_optimize_yaw())
    {
      LOG(INFO)<<"only yaw";
      ceres_pose_control_points.push_back(CeresPose(
          control_points[i], nullptr /* translation_parameterization */,
          nullptr,
          &problem));
      LOG(INFO)<<"knot "<<i<<" before: "<<ceres_pose_control_points[i].ToRigid();
      std::unique_ptr<ceres::LocalParameterization> para(
                      common::make_unique<ceres::AutoDiffLocalParameterization<
                      YawOnlyQuaternionPlus, 4, 1>>());
      problem.AddParameterBlock(ceres_pose_control_points[i].rotation(),4,para.release());
    }
    else
    {
      ceres_pose_control_points.push_back(CeresPose(
          control_points[i], nullptr /* translation_parameterization */,
          nullptr,
          &problem));
      std::unique_ptr<ceres::LocalParameterization> para=
                      common::make_unique<ceres::QuaternionParameterization>();
      problem.AddParameterBlock(ceres_pose_control_points[i].rotation(),4,para.release());
      LOG(INFO)<<"knot "<<i<<" before: "<<ceres_pose_control_points[i].ToRigid();
    }
  }
  std::vector<int> vec = {0,1,2};
  ceres::SubsetParameterization* para_const = new ceres::SubsetParameterization(3,vec);


  /*Old Scan Matcher Modul*/
  problem.AddResidualBlock(
          new ceres:: AutoDiffCostFunction<NurbsControlPointCostFunctor,
          ceres::DYNAMIC, 3,4,3,4,3,4,3,4,3,4>(
              new NurbsControlPointCostFunctor(
                  options_.occupied_space_weight(0)/std::sqrt(static_cast<double>(number_of_residuals)),
                  hybrid_grid[0],
                  begin,
                  end,
                  range_data_vec,
                  ray_tracer_),
                  number_of_residuals),
          nullptr, ceres_pose_control_points[0].translation(), ceres_pose_control_points[0].rotation()
          , ceres_pose_control_points[1].translation(), ceres_pose_control_points[1].rotation()
          , ceres_pose_control_points[2].translation(), ceres_pose_control_points[2].rotation()
          , ceres_pose_control_points[3].translation(), ceres_pose_control_points[3].rotation()
          , ceres_pose_control_points[4].translation(), ceres_pose_control_points[4].rotation());
  problem.AddResidualBlock(
          new ceres:: AutoDiffCostFunction<NurbsControlPointCostFunctor,
          ceres::DYNAMIC, 3,4,3,4,3,4,3,4,3,4>(
              new NurbsControlPointCostFunctor(
                  options_.occupied_space_weight(1)/std::sqrt(static_cast<double>(number_of_residuals)),
                  hybrid_grid[1],
                  begin,
                  end,
                  range_data_vec,
                  ray_tracer_),
                  number_of_residuals),
          nullptr, ceres_pose_control_points[0].translation(), ceres_pose_control_points[0].rotation()
          , ceres_pose_control_points[1].translation(), ceres_pose_control_points[1].rotation()
          , ceres_pose_control_points[2].translation(), ceres_pose_control_points[2].rotation()
          , ceres_pose_control_points[3].translation(), ceres_pose_control_points[3].rotation()
          , ceres_pose_control_points[4].translation(), ceres_pose_control_points[4].rotation());

  /*Linear Acceleration for IMU Cost Functor*/
//  LOG(INFO)<<"residuals:"<<3 * linear_acc_data_vector.size();
//  problem.AddResidualBlock(
//          new ceres:: AutoDiffCostFunction<NurbsTranslationAccelerationDeltaFunctor,
//          ceres::DYNAMIC, 3,4,3,4,3,4,3,4,3,4>(
//              new NurbsTranslationAccelerationDeltaFunctor(
//                  2.0/std::sqrt(static_cast<double>(linear_acc_data_vector.size())),//0.17/std::sqrt(static_cast<double>(3 * linear_acc_data_vector.size())),
//                  linear_acc_data_vector,
//                  begin,
//                  end),
//                  3 * linear_acc_data_vector.size()),
//          nullptr, ceres_pose_control_points[0].translation(), ceres_pose_control_points[0].rotation()
//          , ceres_pose_control_points[1].translation(), ceres_pose_control_points[1].rotation()
//          , ceres_pose_control_points[2].translation(), ceres_pose_control_points[2].rotation()
//          , ceres_pose_control_points[3].translation(), ceres_pose_control_points[3].rotation()
//          , ceres_pose_control_points[4].translation(), ceres_pose_control_points[4].rotation());

  problem.AddResidualBlock(
          new ceres:: AutoDiffCostFunction<NurbsRotationVelocityDeltaFunctor,
          ceres::DYNAMIC, 3,4,3,4,3,4,3,4,3,4>(
              new NurbsRotationVelocityDeltaFunctor(
                  5.0/std::sqrt(static_cast<double>(imu_angular_vel_data_vector.size())),//0.17/std::sqrt(static_cast<double>(3 * linear_acc_data_vector.size())),
                  imu_angular_vel_data_vector,
                  begin,
                  end),
                  3 * imu_angular_vel_data_vector.size()),
          nullptr, ceres_pose_control_points[0].translation(), ceres_pose_control_points[0].rotation()
          , ceres_pose_control_points[1].translation(), ceres_pose_control_points[1].rotation()
          , ceres_pose_control_points[2].translation(), ceres_pose_control_points[2].rotation()
          , ceres_pose_control_points[3].translation(), ceres_pose_control_points[3].rotation()
          , ceres_pose_control_points[4].translation(), ceres_pose_control_points[4].rotation());
  /*new Scan Matcher */
//  freeSpaceEstimator(ceres_pose_control_points,
//                     problem,
//                     number_of_residuals,
//                     hybrid_grid,
//                     decay_grid,
//                     begin,
//                     end,
//                     range_data_vec);

  /*Keep First Control point Constant*/
  problem.SetParameterBlockConstant(ceres_pose_control_points[0].translation());
  problem.SetParameterBlockConstant(ceres_pose_control_points[0].rotation());

  /*Control Point Delta Functor*/
  for (int i = 0; i < control_points.size(); i++){
//    problem.AddResidualBlock(
//        new ceres::AutoDiffCostFunction<TranslationDeltaCostFunctor, 3, 3>(
//            new TranslationDeltaCostFunctor(options_.translation_weight()/std::sqrt(control_points.size()),
//                                            ceres_pose_control_points[i].ToRigid())),
//    nullptr, ceres_pose_control_points[i].translation());
     problem.AddResidualBlock(
         new ceres::AutoDiffCostFunction<RotationDeltaCostFunctor, 3, 4>(
             new RotationDeltaCostFunctor(options_.rotation_weight()/std::sqrt(control_points.size()),
                                          ceres_pose_control_points[i].ToRigid().rotation())),
         nullptr, ceres_pose_control_points[i].rotation());
  }
  LOG(INFO)<<"finished";
  ceres::Solve(ceres_solver_options_, &problem, summary);
  LOG(INFO)<<"solve";
  LOG(INFO)<<summary->FullReport();
  for (int i = 0; i < ceres_pose_control_points.size(); i++){
    control_points[i] = ceres_pose_control_points[i].ToRigid();
    Eigen::Quaterniond quat_norm = control_points[i].rotation();
    quat_norm.normalize();
    Eigen::Matrix3d rot_mat = Eigen::Matrix3d(quat_norm);
    Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
    LOG(INFO)<<" roll: "<<angles(0)*180/M_PI<<" pitch: "<<angles(1)*180/M_PI<<" yaw: "<<angles(2)*180/M_PI;
    control_points[i] = transform::Rigid3d(control_points[i].translation(), quat_norm);
    LOG(INFO)<<"knot "<<i<<" after: "<<ceres_pose_control_points[i].ToRigid();
  }
  LOG(INFO)<<"Time Batch Duration:"<<common::ToSeconds(end-begin);

}
void CeresScanMatcher::freeSpaceEstimator(std::vector<CeresPose>& ceres_pose_control_points,
                                          ceres::Problem& problem,
                                          const int& number_of_residuals,
                                          const std::vector<const HybridGrid*>& hybrid_grid,
                                          const std::vector<const HybridDecayGrid*>& decay_grid,
                                          const common::Time& begin,
                                          const common::Time& end,
                                          const std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec)
{
  for (int i = 0; i < hybrid_grid.size() - 1; i++)
  {
    problem.AddResidualBlock(
              new ceres::AutoDiffCostFunction<NurbsControlPointFreeSpaceCostFunctor, ceres::DYNAMIC,3,4,3,4,3,4,3,4,3,4>(
                  new NurbsControlPointFreeSpaceCostFunctor(
                      options_.occupied_space_weight(i)/std::sqrt(static_cast<double>(number_of_residuals)),
                      *hybrid_grid[i],
                      decay_grid[i],
                      begin,
                      end,
                      range_data_vec,
                      ray_tracer_,
                      options_.ray_tracer_line_size()), number_of_residuals),
              nullptr, ceres_pose_control_points[0].translation(), ceres_pose_control_points[0].rotation()
              , ceres_pose_control_points[1].translation(), ceres_pose_control_points[1].rotation()
              , ceres_pose_control_points[2].translation(), ceres_pose_control_points[2].rotation()
              , ceres_pose_control_points[3].translation(), ceres_pose_control_points[3].rotation()
              , ceres_pose_control_points[4].translation(), ceres_pose_control_points[4].rotation());
    }
}

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer
