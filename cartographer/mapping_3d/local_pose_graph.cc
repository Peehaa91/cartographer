/*
 * local_pose_graph.cc
 *
 *  Created on: Jul 15, 2017
 *      Author: schnattinger
 */
#include "cartographer/mapping_3d/local_pose_graph.h"
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace cartographer {
namespace mapping_3d {



proto::LocalPoseGraphOptions CreateLocalPoseGraphOptions(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::LocalPoseGraphOptions options;
  options.set_knot_type(parameter_dictionary->GetString("knot_type"));
  options.set_weight_type(parameter_dictionary->GetString("weight_type"));
  options.set_debug(parameter_dictionary->GetBool("debug"));
  return options;
}

LocalPoseGraph::LocalPoseGraph(const proto::LocalPoseGraphOptions& options)
    : options_(options),
      input_size_(1),
      output_size_(1),
      knot_type_(KnotType::UNIFORM),
      weight_type_(WeightType::NON_RATIONAL) {

}

void LocalPoseGraph::createPoseAndScanFromSpline(std::vector<PoseEstimate>& control_point_vec,
                                                std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
                                                std::vector<PoseAndRangeData>& pose_and_cloud_vec)
{
  boost::multi_array<std::array<double, output_dim>, input_dim> points(
      boost::extents[control_point_vec.size()]);
  for (int i = 0; i < control_point_vec.size(); i++) {
    points[i][0] = control_point_vec[i].pose.translation().x();
    points[i][1] = control_point_vec[i].pose.translation().y();
    points[i][2] = control_point_vec[i].pose.translation().z();
    Eigen::Quaterniond quat = control_point_vec[i].pose.rotation();
    points[i][3] = control_point_vec[i].pose.rotation().x();
    points[i][4] = control_point_vec[i].pose.rotation().y();
    points[i][5] = control_point_vec[i].pose.rotation().z();
    points[i][6] = control_point_vec[i].pose.rotation().w();
  }
  NurbsKnots<double, input_dim, output_dim, knot_type> knots;
  knots.create(degree, points.shape());

  //create weight vector
  NurbsWeights<double, input_dim, output_dim, weight_type> weights;
  weights.create(points.shape());

  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<
      double, input_dim, output_dim, knot_type, weight_type>();
  nurbs.init(degree, knots, weights, points);
  std::array<std::vector<double>, input_dim> knot_vec = knots.getKnot();
  for (int i = 0; i < knot_vec[0].size(); i++)
    LOG(INFO)<<"knot: "<<knot_vec[0][i];

  points= nurbs.getPoints();
  common::Time begin = control_point_vec.front().time;
  common::Time end = control_point_vec.back().time;

  int numSamplePoints = range_data_vec.size();
  for (int i = 0; i < numSamplePoints; ++i) {
    double point = common::ToSeconds(range_data_vec[i].first - begin)
        / common::ToSeconds(end - begin);
    std::array<double, output_dim> output;
    nurbs.getPoint(&point, output.begin());
    transform::Rigid3d pose = transform::Rigid3d(
        Eigen::Vector3d(output[0], output[1], output[2]),
        Eigen::Quaterniond(output[6], output[3], output[4], output[5]));
    Eigen::Vector3d euler = Eigen::Quaterniond(output[6], output[3], output[4], output[5]).toRotationMatrix().eulerAngles(0, 1, 2);
    pose_and_cloud_vec.push_back(PoseAndRangeData(range_data_vec[i].first, pose, range_data_vec[i].second));
  }
}
std::vector<PoseAndRangeData> LocalPoseGraph::createSplineFromControlVector(
    std::vector<PoseEstimate>& control_point_vec,
    std::vector<std::pair<common::Time, sensor::RangeData>>& range_data_vec,
    std::vector<std::pair<common::Time, Eigen::Vector3d>>& imu_angular_vel_data_vector,
    std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector,
    const HybridGrid* hybrid_grid,
    const std::unique_ptr<scan_matching::CeresScanMatcher>& ceres_scan_matcher,
    const HybridDecayGrid* decay_grid) {
  if (options_.debug())
  {
    std::vector<PoseAndRangeData> init_pose_and_cloud_vec;
    createPoseAndScanFromSpline(control_point_vec, range_data_vec,
                                init_pose_and_cloud_vec);
    writeSplineInFile(control_point_vec, init_pose_and_cloud_vec, std::to_string(nurbs_number) + "_init");
  }
  Eigen::Matrix3d rot_mat = Eigen::Matrix3d(control_point_vec[0].pose.rotation());
  Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
//  LOG(INFO)<<"first: "<<pose_vec[0].pose.translation()<<angles * 180/M_PI;
  rot_mat = Eigen::Matrix3d(control_point_vec[control_point_vec.size() - 1].pose.rotation());
  angles = rot_mat.eulerAngles(0, 1, 2);
//  LOG(INFO)<<"last: "<<pose_vec[pose_vec.size()-1].pose<<angles * 180/M_PI;
  std::vector<transform::Rigid3d> knot_vector;
  boost::multi_array<std::array<double, output_dim>, input_dim> points(
      boost::extents[control_point_vec.size()]);
  for (int i = 0; i < control_point_vec.size(); i++) {
//    LOG(INFO)<<"pose "<<i<<": "<<pose_vec[i].pose;
    points[i][0] = control_point_vec[i].pose.translation().x();
    points[i][1] = control_point_vec[i].pose.translation().y();
    points[i][2] = control_point_vec[i].pose.translation().z();
    Eigen::Quaterniond quat = control_point_vec[i].pose.rotation();
    points[i][3] = control_point_vec[i].pose.rotation().x();
    points[i][4] = control_point_vec[i].pose.rotation().y();
    points[i][5] = control_point_vec[i].pose.rotation().z();
    points[i][6] = control_point_vec[i].pose.rotation().w();
    knot_vector.push_back(control_point_vec[i].pose);
//    rot_mat = Eigen::Matrix3d(quat);
//    angles = rot_mat.eulerAngles(0, 1, 2);
//    if (i > 0 && (fabs(angles(0) - points[i-1][3]) > 2))
//      LOG(INFO)<<"orientation before: "<<pose_vec[i-1].pose<<" orientation now: "<<pose_vec[i].pose;
//    points[i][3] = angles(0);
//    points[i][4] = angles(1);
//    points[i][5] = angles(2);

  }
//  LOG(INFO)<<"pose before: "<<pose_vec.back().pose;
//  LOG(INFO)<<"size: "<<range_data_vec.size();
//  NurbsKnots<double, input_dim, output_dim, knot_type> knots;
//  knots.create(degree, points.shape());
//
//  //create weight vector
//  NurbsWeights<double, input_dim, output_dim, weight_type> weights;
//  weights.create(points.shape());
//
//  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<
//      double, input_dim, output_dim, knot_type, weight_type>();
//  nurbs.init(degree, knots, weights, points);

  common::Time begin = control_point_vec.front().time;
  common::Time end = control_point_vec.back().time;

  ceres::Solver::Summary summary;

//  ceres_scan_matcher->Match(hybrid_grid,
//                            begin,
//                            end,
//                            range_data_vec,
//                            nurbs,
//                            &nurbs,
//                            &summary);
  std::vector<double> weight_vec = {1, 1, 1, 1, 1};
  ceres_scan_matcher->Match(hybrid_grid,
                            decay_grid,
                            begin,
                            end,
                            range_data_vec,
                            imu_angular_vel_data_vector,
                            linear_acc_data_vector,
                            knot_vector,
                            weight_vec,
                            &summary);


  for (int i = 0; i < control_point_vec.size(); i++) {
//    LOG(INFO)<<"pose "<<i<<": "<<pose_vec[i].pose;
    points[i][0] = knot_vector[i].translation().x();
    points[i][1] = knot_vector[i].translation().y();
    points[i][2] = knot_vector[i].translation().z();
    Eigen::Quaterniond quat = control_point_vec[i].pose.rotation();
    points[i][3] = knot_vector[i].rotation().x();
    points[i][4] = knot_vector[i].rotation().y();
    points[i][5] = knot_vector[i].rotation().z();
    points[i][6] = knot_vector[i].rotation().w();
//    rot_mat = Eigen::Matrix3d(quat);
//    angles = rot_mat.eulerAngles(0, 1, 2);
//    if (i > 0 && (fabs(angles(0) - points[i-1][3]) > 2))
//      LOG(INFO)<<"orientation before: "<<pose_vec[i-1].pose<<" orientation now: "<<pose_vec[i].pose;
//    points[i][3] = angles(0);
//    points[i][4] = angles(1);
//    points[i][5] = angles(2);

  }
  std::vector<PoseEstimate> control_points;
  for (int i = 0; i < control_point_vec.size(); i++)
  {
    transform::Rigid3d pose = transform::Rigid3d(
        Eigen::Vector3d(points[i][0], points[i][1], points[i][2]),
            Eigen::Quaterniond(points[i][6], points[i][3], points[i][4], points[i][5]));
    PoseEstimate pose_est;
    pose_est.time = control_point_vec[i].time;
    pose_est.pose = pose;
    control_points.push_back(pose_est);
  }
  std::vector<PoseAndRangeData> pose_and_cloud_vec;
  createPoseAndScanFromSpline(control_points,
                              range_data_vec,
                              pose_and_cloud_vec);
  writeSplineInFile(control_points ,pose_and_cloud_vec, std::to_string(nurbs_number));
  nurbs_number++;
  return pose_and_cloud_vec;

}

void LocalPoseGraph::writeSplineInFile(std::vector<PoseEstimate>& control_point_vec,
                                       std::vector<PoseAndRangeData> & range_data_vec,
                                       std::string file_name)
{
  struct stat myStat;
  const char* path = "/home/schnattinger/.ros/nurbs";
  int result = stat(path, &myStat);
  LOG(INFO)<<"output of stat: "<<result;
  if (result == -1){
    LOG(INFO)<<"created directory: "<<path;
    mkdir(path, 0700);
  }
  std::string file_path = std::string(path) + std::string("/") + file_name + std::string(".txt");
  std::ofstream nurbs_file(file_path.c_str());
  nurbs_file<<"knots:"<<std::endl;
  for (PoseEstimate& pose : control_point_vec)
  {
    nurbs_file<<"time: "<<pose.time<<" "<<pose.pose<<std::endl;
  }
  nurbs_file<<"points:"<<std::endl;
  for (PoseAndRangeData pose_and_range : range_data_vec)
  {
    nurbs_file<<"time: "<<pose_and_range.time<<" "<<pose_and_range.pose<<std::endl;
  }
  nurbs_file.close();


}
}
}
