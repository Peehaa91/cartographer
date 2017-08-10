/*
 * local_pose_graph.cc
 *
 *  Created on: Jul 15, 2017
 *      Author: schnattinger
 */
#include "cartographer/mapping_3d/local_pose_graph.h"

namespace cartographer {
namespace mapping_3d {

proto::LocalPoseGraphOptions CreateLocalPoseGraphOptions(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::LocalPoseGraphOptions options;
  options.set_knot_type(parameter_dictionary->GetString("knot_type"));
  options.set_weight_type(parameter_dictionary->GetString("weight_type"));
  return options;
}

LocalPoseGraph::LocalPoseGraph(const proto::LocalPoseGraphOptions& options)
    : options_(options),
      input_size_(1),
      output_size_(1),
      knot_type_(KnotType::UNIFORM),
      weight_type_(WeightType::NON_RATIONAL) {

}

std::vector<PoseAndRangeData> LocalPoseGraph::createSplineFromControlVector(
    std::vector<PoseEstimate>& pose_vec,
    std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
    const HybridGrid* hybrid_grid,
    const std::unique_ptr<scan_matching::CeresScanMatcher>& ceres_scan_matcher) {
  const int input_dim = 1;
  const int output_dim = 6;
  const int degree = 2;
  const KnotType knot_type = KnotType::UNIFORM;
  const WeightType weight_type = WeightType::RATIONAL;
  Eigen::Matrix3d rot_mat = Eigen::Matrix3d(pose_vec[0].pose.rotation());
  Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
//  LOG(INFO)<<"first: "<<pose_vec[0].pose.translation()<<angles * 180/M_PI;
  rot_mat = Eigen::Matrix3d(pose_vec[pose_vec.size() - 1].pose.rotation());
  angles = rot_mat.eulerAngles(0, 1, 2);
//  LOG(INFO)<<"last: "<<pose_vec[pose_vec.size()-1].pose<<angles * 180/M_PI;
  boost::multi_array<std::array<double, output_dim>, input_dim> points(
      boost::extents[pose_vec.size()]);
  for (int i = 0; i < pose_vec.size(); i++) {
    points[i][0] = pose_vec[i].pose.translation().x();
    points[i][1] = pose_vec[i].pose.translation().y();
    points[i][2] = pose_vec[i].pose.translation().z();
    Eigen::Quaterniond quat = pose_vec[i].pose.rotation();
    rot_mat = Eigen::Matrix3d(quat);
    angles = rot_mat.eulerAngles(0, 1, 2);
    points[i][3] = angles(0);
    points[i][4] = angles(1);
    points[i][5] = angles(2);

  }

  LOG(INFO)<<"size: "<<range_data_vec.size();
  NurbsKnots<double, input_dim, output_dim, knot_type> knots;
  knots.create(degree, points.shape());

  //create weight vector
  NurbsWeights<double, input_dim, output_dim, weight_type> weights;
  weights.create(points.shape());

  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<
      double, input_dim, output_dim, knot_type, weight_type>();
  nurbs.init(degree, knots, weights, points);


  double minU = nurbs.getMinU(0);
  double maxU = nurbs.getMaxU(0);
  common::Time begin = pose_vec.front().time;
  common::Time end = pose_vec.back().time;

  ceres::Solver::Summary summary;

  ceres_scan_matcher->Match(hybrid_grid,
                            begin,
                            end,
                            range_data_vec,
                            nurbs,
                            &nurbs,
                            &summary);
  boost::multi_array<double, 1>& weights_estimate = nurbs.getWeights().getWeights();
  for (int i = 0; i < weights_estimate.size(); i++) {
    LOG(INFO)<<"weight: "<<weights_estimate[i];
    }
  int numSamplePoints = range_data_vec.size();
  std::vector<PoseAndRangeData> pose_and_cloud_vec;
  for (int i = 0; i < numSamplePoints; ++i) {
//    double point = minU
//        + (maxU - minU) / (double) (numSamplePoints - 1) * (double) i;
    double point = common::ToSeconds(range_data_vec[i].first - begin)
        / common::ToSeconds(end - begin);
    std::array<double, output_dim> output;
    nurbs.getPoint(&point, output.begin());
    transform::Rigid3d pose = transform::Rigid3d(
        Eigen::Vector3d(output[0], output[1], output[2]),
        transform::AngleAxisVectorToRotationQuaternion(
            Eigen::Vector3d(output[3], output[4], output[5])));
    pose_and_cloud_vec.push_back(PoseAndRangeData(range_data_vec[i].first, pose, range_data_vec[i].second));
//    LOG(INFO)<<"u: "<<point<<" x: "<<output[0]<<" y: "<<output[1]<<" z: "<<output[2]<<" roll: "<<output[3]*180/M_PI<<" pitch: "<<output[4]*180/M_PI<<" yaw: "<<output[5]*180/M_PI;
  }
  return pose_and_cloud_vec;

}
}
}
