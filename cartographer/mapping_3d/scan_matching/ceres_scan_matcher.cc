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
#include "cartographer/mapping_3d/scan_matching/nurbs_knot_cost_functor.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/transform/transform.h"
#include "ceres/ceres.h"
#include "glog/logging.h"

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
  for (size_t i = 0; i != point_clouds_and_hybrid_grids.size(); ++i) {
    CHECK_GT(options_.occupied_space_weight(i), 0.);
    const sensor::PointCloud& point_cloud =
        *point_clouds_and_hybrid_grids[i].first;
    const HybridGrid& hybrid_grid = *point_clouds_and_hybrid_grids[i].second;
//    LOG(INFO)<<"ceres before: "<<ceres_pose.ToRigid();
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
//  LOG(INFO)<<"ceres after: "<<ceres_pose.ToRigid();
  sum += summary->total_time_in_seconds;
  counter++;
  LOG(INFO)<<"mean secs: "<<sum/counter;
//  LOG(INFO)<<summary->FullReport();
  *pose_estimate = ceres_pose.ToRigid();
}

void CeresScanMatcher::Match(const HybridGrid* hybrid_grid,
                             const common::Time& begin,
                             const common::Time& end,
                             const std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
                             Nurbs<double, 1, 6, KnotType::UNIFORM, WeightType::RATIONAL>& nurbs,
                             Nurbs<double, 1, 6, KnotType::UNIFORM, WeightType::RATIONAL>* nurbs_estimate,
                             ceres::Solver::Summary* summary)
{
  ceres::Problem problem;
  int number_of_residuals = 0;
  for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec)
    number_of_residuals += data.second.returns.size();
  NurbsWeights<double, 1, 6, WeightType::RATIONAL>& weights = nurbs.getWeights();
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
  boost::multi_array<std::array<double, 6>, 1>& points = nurbs.getPoints();
  problem.AddResidualBlock(
          new ceres:: NumericDiffCostFunction<NurbsKnotCostFunctor,
          ceres::FORWARD,ceres::DYNAMIC, 6, 6, 6, 6, 6>(
              new NurbsKnotCostFunctor(
                  options_.occupied_space_weight(0)/std::sqrt(static_cast<double>(number_of_residuals)),
                  hybrid_grid_copy,
                  begin,
                  end,
                  range_data_vec,
                  &nurbs), ceres::TAKE_OWNERSHIP,
                  number_of_residuals),
          nullptr, points[0].data(), points[1].data(), points[2].data(), points[3].data(), points[4].data());
  for (int i = 0; i < 5; i++)
  {
    problem.AddResidualBlock(
        new ceres::AutoDiffCostFunction<KnotTranslationDeltaCostFunctor, 3, 6>(
            new KnotTranslationDeltaCostFunctor(options_.translation_weight(),
                                            points[i])),
        nullptr, points[i].data());
  }
//  problem.AddResidualBlock(
//          new ceres:: AutoDiffCostFunction<NurbsWeightCostFunctor,
//          ceres::DYNAMIC, 5>(
//              new NurbsWeightCostFunctor(
//                  options_.occupied_space_weight(0)/std::sqrt(static_cast<double>(number_of_residuals)),
//                  hybrid_grid_copy,
//                  begin,
//                  end,
//                  range_data_vec,
//                  nurbs),
//                  number_of_residuals),
//          nullptr, initial_weights.data());
//  for (int i = 0; i < weights_vec.size(); i++) {
//    problem.SetParameterLowerBound(initial_weights.data(), i, 0);
//    problem.SetParameterUpperBound(initial_weights.data(), i, 1);

//  }
  std::unique_ptr<ceres::LocalParameterization> parametrization = common::make_unique<ceres:: HomogeneousVectorParameterization>(5);
//  problem.SetParameterization(initial_weights.data(),parametrization.release());
  LOG(INFO)<<"number of threads:"<<ceres_solver_options_.num_threads;
  ceres::Solve(ceres_solver_options_, &problem, summary);
  LOG(INFO)<<summary->FullReport();
  boost::multi_array<double, 1>& weights_estimate = nurbs_estimate->getWeights().getWeights();
  for (int i = 0; i < initial_weights.size(); i++) {
    weights_estimate[i] = initial_weights[i];
    }

}

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer
