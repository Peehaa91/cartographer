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
    const std::vector<const HybridGrid*>& hybrid_grid,
    const std::unique_ptr<scan_matching::CeresScanMatcher>& ceres_scan_matcher,
    const std::vector<const HybridDecayGrid*>& decay_grid) {
  if (options_.debug())
  {
    std::vector<PoseAndRangeData> init_pose_and_cloud_vec;
    createPoseAndScanFromSpline(control_point_vec, range_data_vec,
                                init_pose_and_cloud_vec);
//    writeSplineInFile(control_point_vec, init_pose_and_cloud_vec, std::to_string(nurbs_number) + "_init");
  }
  Eigen::Matrix3d rot_mat = Eigen::Matrix3d(control_point_vec[0].pose.rotation());
  Eigen::Vector3d angles = rot_mat.eulerAngles(0, 1, 2);
//  LOG(INFO)<<"first: "<<pose_vec[0].pose.translation()<<angles * 180/M_PI;
  rot_mat = Eigen::Matrix3d(control_point_vec[control_point_vec.size() - 1].pose.rotation());
  angles = rot_mat.eulerAngles(0, 1, 2);
//  LOG(INFO)<<"last: "<<pose_vec[pose_vec.size()-1].pose<<angles * 180/M_PI;
  std::vector<transform::Rigid3d> control_points_pose_vec;
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
    control_points_pose_vec.push_back(control_point_vec[i].pose);

  }


  common::Time begin = control_point_vec.front().time;
  common::Time end = control_point_vec.back().time;

  ceres::Solver::Summary summary;

  std::vector<double> weight_vec = {1, 1, 1, 1, 1};
  ceres_scan_matcher->MatchSplineWithFreeSpace(hybrid_grid,
                            decay_grid,
                            begin,
                            end,
                            range_data_vec,
                            imu_angular_vel_data_vector,
                            linear_acc_data_vector,
                            control_points_pose_vec,
                            weight_vec,
                            &summary);


  for (int i = 0; i < control_point_vec.size(); i++) {
//    LOG(INFO)<<"pose "<<i<<": "<<pose_vec[i].pose;
    points[i][0] = control_points_pose_vec[i].translation().x();
    points[i][1] = control_points_pose_vec[i].translation().y();
    points[i][2] = control_points_pose_vec[i].translation().z();
    Eigen::Quaterniond quat = control_point_vec[i].pose.rotation();
    points[i][3] = control_points_pose_vec[i].rotation().x();
    points[i][4] = control_points_pose_vec[i].rotation().y();
    points[i][5] = control_points_pose_vec[i].rotation().z();
    points[i][6] = control_points_pose_vec[i].rotation().w();

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
  writeSplineInFile(control_points ,pose_and_cloud_vec, std::to_string(nurbs_number), linear_acc_data_vector);
  control_point_vec = control_points;
  nurbs_number++;
  return pose_and_cloud_vec;

}
void LocalPoseGraph::createDerivativeSplines(std::vector<PoseEstimate>& control_point_vec,
                                             std::vector<PoseAndRangeData> & range_data_vec,
                                             std::vector<std::pair<common::Time, transform::Rigid3d>>& velocity_data,
                                             std::vector<std::pair<common::Time, transform::Rigid3d>>& acceleration_data,
                                             std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector,
                                             std::vector<std::pair<common::Time, Eigen::Vector3d>>& imu_res_vec)
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
  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurbs_first_deriv = Nurbs<
      double, input_dim, output_dim, knot_type, weight_type>();
  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurbs_second_deriv = Nurbs<
      double, input_dim, output_dim, knot_type, weight_type>();
  NurbsKnots<double, input_dim, output_dim, knot_type> knots_first_deriv = NurbsKnots<
      double, input_dim, output_dim, knot_type>();
  NurbsKnots<double, input_dim, output_dim, knot_type> knots_second_deriv = NurbsKnots<
      double, input_dim, output_dim, knot_type>();
  NurbsWeights<double, input_dim, output_dim, weight_type> weights_first_deriv =
      NurbsWeights<double, input_dim, output_dim, weight_type>();
  NurbsWeights<double, input_dim, output_dim, weight_type> weights_second_deriv =
      NurbsWeights<double, input_dim, output_dim, weight_type>();

  std::array<std::vector<double>, input_dim> knot_vec = knots.getKnot();
  boost::multi_array<std::array<double, output_dim>, input_dim> points_first_deriv(
          boost::extents[control_point_vec.size() -1]);
      for (int i = 0; i < control_point_vec.size() -1; ++i) {
        points_first_deriv[i][0] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][0] - points[i][0]);
        points_first_deriv[i][1] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][1] - points[i][1]);
        points_first_deriv[i][2] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][2] - points[i][2]);
        points_first_deriv[i][3] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][3] - points[i][3]);
        points_first_deriv[i][4] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][4] - points[i][4]);
        points_first_deriv[i][5] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][5] - points[i][5]);
        points_first_deriv[i][6] = (degree)/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][6] - points[i][6]);
      }
      knots_first_deriv.create(degree - 1, points_first_deriv.shape());
      weights_first_deriv.create(points_first_deriv.shape());
      nurbs_first_deriv.init(degree - 1, knots_first_deriv, weights_first_deriv, points_first_deriv);
      common::Time begin = control_point_vec.front().time;
      common::Time end = control_point_vec.back().time;
      int numSamplePoints = range_data_vec.size();
      for (int i = 0; i < numSamplePoints; ++i) {
        double point = common::ToSeconds(range_data_vec[i].time - begin)
            / common::ToSeconds(end - begin);
        std::array<double, output_dim> output;
        nurbs_first_deriv.getPoint(&point, output.begin());
        transform::Rigid3d pose = transform::Rigid3d(
            Eigen::Vector3d(output[0], output[1], output[2]),
            Eigen::Quaterniond(output[6], output[3], output[4], output[5]));
        velocity_data.push_back(std::make_pair(range_data_vec[i].time, pose));
      }



      std::array<std::vector<double>, input_dim> knot_first_vec = knots_first_deriv.getKnot();

      boost::multi_array<std::array<double, output_dim>, input_dim> points_second_deriv(
          boost::extents[control_point_vec.size() -2]);
      for (int i = 0; i < control_point_vec.size() -2; ++i) {
        points_second_deriv[i][0] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][0] - points_first_deriv[i][0]);
        points_second_deriv[i][1] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][1] - points_first_deriv[i][1]);
        points_second_deriv[i][2] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][2] - points_first_deriv[i][2]);
        points_second_deriv[i][3] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][3] - points_first_deriv[i][3]);
        points_second_deriv[i][4] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][4] - points_first_deriv[i][4]);
        points_second_deriv[i][5] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][5] - points_first_deriv[i][5]);
        points_second_deriv[i][6] = (degree -1)/(knot_first_vec[0][i + degree -1 + 1] - knot_first_vec[0][i + 1])
            * (points_first_deriv[i + 1][6] - points_first_deriv[i][6]);
      }
      weights_second_deriv.create(points_second_deriv.shape());
      knots_second_deriv.create(degree - 2, points_second_deriv.shape());
      nurbs_second_deriv.init(degree - 2, knots_second_deriv, weights_second_deriv, points_second_deriv);
      for (int i = 0; i < numSamplePoints; ++i) {
        double point = common::ToSeconds(range_data_vec[i].time - begin)
            / common::ToSeconds(end - begin);
        std::array<double, output_dim> output;
        nurbs_second_deriv.getPoint(&point, output.begin());
        transform::Rigid3d pose = transform::Rigid3d(
            Eigen::Vector3d(output[0], output[1], output[2]),
            Eigen::Quaterniond(output[6], output[3], output[4], output[5]));
        acceleration_data.push_back(std::make_pair(range_data_vec[i].time, pose));
      }


      for (int i = 0; i < linear_acc_data_vector.size(); ++i) {
        double point = common::ToSeconds(linear_acc_data_vector[i].first - begin)
            / common::ToSeconds(end - begin);
        std::array<double, output_dim> output;
        nurbs.getPoint(&point, output.begin());
        transform::Rigid3d pose = transform::Rigid3d(
            Eigen::Vector3d(output[0], output[1], output[2]),
            Eigen::Quaterniond(output[6], output[3], output[4], output[5]));
        Eigen::Vector3d gravity(0, 0, 9.81);
        Eigen::Vector3d gravity_rotated = pose.rotation().inverse() * gravity;
        Eigen::Vector3d acceleration = linear_acc_data_vector[i].second - gravity_rotated;
        imu_res_vec.push_back(std::make_pair(linear_acc_data_vector[i].first, acceleration));
      }

}

void LocalPoseGraph::writeSplineInFile(std::vector<PoseEstimate>& control_point_vec,
                                       std::vector<PoseAndRangeData> & range_data_vec,
                                       std::string file_name,
                                       std::vector<std::pair<common::Time, Eigen::Vector3d>>& linear_acc_data_vector)
{
  std::vector<std::pair<common::Time, transform::Rigid3d>> velocity_data;
  std::vector<std::pair<common::Time, transform::Rigid3d>> acceleration_data;
  std::vector<std::pair<common::Time, Eigen::Vector3d>> imu_res_vec;
  createDerivativeSplines(control_point_vec,
                           range_data_vec,
                           velocity_data,
                           acceleration_data,
                           linear_acc_data_vector,
                           imu_res_vec);
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
  file_path = std::string(path) + std::string("/") + file_name + std::string("_trans.txt");
  std::ofstream nurbs_trans_file(file_path.c_str());
  for (PoseAndRangeData pose_and_range : range_data_vec)
  {
    nurbs_trans_file<<"time: "<< common::ToSeconds(pose_and_range.time - control_point_vec[0].time)<<" x: "<<pose_and_range.pose.translation().x()
    <<" y: "<<pose_and_range.pose.translation().y()
    <<" z: "<<pose_and_range.pose.translation().z()<<std::endl;
  }
  nurbs_trans_file.close();

  file_path = std::string(path) + std::string("/") + file_name + std::string("_vel.txt");
  std::ofstream nurbs_vel_file(file_path.c_str());
  for (std::pair<common::Time, transform::Rigid3d> vel_data : velocity_data)
  {
    nurbs_vel_file<<"time: "<< common::ToSeconds(vel_data.first- control_point_vec[0].time)<<" x: "<<vel_data.second.translation().x()
    <<" y: "<<vel_data.second.translation().y()
    <<" z: "<<vel_data.second.translation().z()<<std::endl;
  }
  nurbs_vel_file.close();

  file_path = std::string(path) + std::string("/") + file_name + std::string("_acc.txt");
  std::ofstream nurbs_acc_file(file_path.c_str());
  for (std::pair<common::Time, transform::Rigid3d> acc_data : acceleration_data)
  {
    nurbs_acc_file<<"time: "<< common::ToSeconds(acc_data.first- control_point_vec[0].time)<<" x: "<<acc_data.second.translation().x()
    <<" y: "<<acc_data.second.translation().y()
    <<" z: "<<acc_data.second.translation().z()<<std::endl;
  }
  nurbs_acc_file.close();

  file_path = std::string(path) + std::string("/") + file_name + std::string("_acc_imu.txt");
  std::ofstream nurbs_acc_imu_file(file_path.c_str());
  for (std::pair<common::Time, Eigen::Vector3d> acc_data : linear_acc_data_vector)
  {
    nurbs_acc_imu_file<<"time: "<< common::ToSeconds(acc_data.first- control_point_vec[0].time)<<" x: "<<acc_data.second.x()
    <<" y: "<<acc_data.second.y()
    <<" z: "<<acc_data.second.z()<<std::endl;
  }
  nurbs_acc_imu_file.close();
  file_path = std::string(path) + std::string("/") + file_name + std::string("res_acc_imu.txt");
  std::ofstream nurbs_res_acc_imu_file(file_path.c_str());
  for (std::pair<common::Time, Eigen::Vector3d> acc_data : imu_res_vec)
  {
    nurbs_res_acc_imu_file<<"time: "<< common::ToSeconds(acc_data.first- control_point_vec[0].time)<<" x: "<<acc_data.second.x()
    <<" y: "<<acc_data.second.y()
    <<" z: "<<acc_data.second.z()<<std::endl;
  }
  nurbs_res_acc_imu_file.close();




}
}
}
