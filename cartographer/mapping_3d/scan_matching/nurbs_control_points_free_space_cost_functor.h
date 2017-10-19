/*
 * nurbs_free_space_cost_function.h
 *
 *  Created on: Sep 8, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_CONTROL_POINTS_FREE_SPACE_COST_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_CONTROL_POINTS_FREE_SPACE_COST_FUNCTOR_H_

#include "Eigen/Core"
#include "cartographer/mapping_3d/hybrid_grid.h"
#include "cartographer/mapping_3d/scan_matching/interpolated_decay_grid.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/transform/rigid_transform.h"
#include "cartographer/transform/transform.h"
#include "cartographer/mapping_3d/scan_matching/nurbs.h"
#include "cartographer/mapping_3d/ray_tracer.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {

// Computes the cost of inserting occupied space described by the point cloud
// into the map. The cost increases with the amount of free space that would be
// replaced by occupied space.
class NurbsControlPointFreeSpaceCostFunctor {
 public:
  // Creates an OccupiedSpaceCostFunctor using the specified grid, 'rotation' to
  // add to all poses, and point cloud.
  NurbsControlPointFreeSpaceCostFunctor(
      const double scaling_factor,
      const HybridGrid& hybrid_grid,
      const HybridDecayGrid* decay_grid,
      const common::Time& begin,
      const common::Time& end,
      const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec,
      const std::shared_ptr<RayTracer>& ray_tracer, int line_size)
      : scaling_factor_(scaling_factor),
        range_data_vec_(range_data_vec),
        interpolated_decay_grid_(*decay_grid),
        interpolated_grid_(hybrid_grid),
        hybrid_grid_(&hybrid_grid),
        decay_grid_(decay_grid),
        begin_(begin),
        end_(end),
        ray_tracer_(ray_tracer),
        line_size_(line_size) {
  }

  NurbsControlPointFreeSpaceCostFunctor(
      const NurbsControlPointFreeSpaceCostFunctor&) = delete;
  NurbsControlPointFreeSpaceCostFunctor& operator=(
      const NurbsControlPointFreeSpaceCostFunctor&) = delete;

  template<typename T>
  double getScalar(const T& jet) const {
    return ceres::Jet<double, 35>(jet).a;
  }

  template<typename T>
  Eigen::Matrix<double, 35, 1> getCov(const T jet) const {
    return ceres::Jet<double, 35>(jet).v;
  }

  template<typename T>
  bool operator()(const T* const point_1_trans, const T* const point_1_rotation,
                  const T* const point_2_trans, const T* const point_2_rotation,
                  const T* const point_3_trans, const T* const point_3_rotation,
                  const T* const point_4_trans, const T* const point_4_rotation,
                  const T* const point_5_trans, const T* const point_5_rotation,
                  T* const residual) const {

    std::vector<const T*> vec_trans = { point_1_trans, point_2_trans,
        point_3_trans, point_4_trans, point_5_trans };
    std::vector<const T*> vec_rot = { point_1_rotation, point_2_rotation,
        point_3_rotation, point_4_rotation, point_5_rotation };
    Nurbs<T, input_dim, output_dim, knot_type, weight_type> nurbs = Nurbs<T,
        input_dim, output_dim, knot_type, weight_type>();
    NurbsKnots<T, input_dim, output_dim, knot_type> knots = NurbsKnots<T,
        input_dim, output_dim, knot_type>();
    NurbsWeights<T, input_dim, output_dim, weight_type> weights = NurbsWeights<
        T, input_dim, output_dim, weight_type>();
    boost::multi_array<std::array<T, output_dim>, input_dim> points(
        boost::extents[vec_trans.size()]);
    for (int i = 0; i < vec_trans.size(); i++) {
      points[i][0] = vec_trans[i][0];
      points[i][1] = vec_trans[i][1];
      points[i][2] = vec_trans[i][2];
      points[i][3] = vec_rot[i][1];
      points[i][4] = vec_rot[i][2];
      points[i][5] = vec_rot[i][3];
      points[i][6] = vec_rot[i][0];

    }
    knots.create(degree, points.shape());
    weights.create(points.shape());
    nurbs.init(degree, knots, weights, points);
    points = nurbs.getPoints();
    int numSamplePoints = range_data_vec_.size();
    sensor::PointCloud complete_pointcloud;
    int counter = 0;
    for (int i = 0; i < numSamplePoints; ++i) {
      //    double point = minU
      //        + (maxU - minU) / (double) (numSamplePoints - 1) * (double) i;
      double point = ((common::ToSeconds(range_data_vec_[i].first - begin_)
          / common::ToSeconds(end_ - begin_)));
      std::array<T, output_dim> output;
      nurbs.getPoint(&point, output.begin());
      const transform::Rigid3<T> transform(
          Eigen::Matrix<T, 3, 1>(output[0], output[1], output[2]),
          Eigen::Quaternion<T>(output[6], output[3], output[4], output[5]));
      for (Eigen::Vector3f point : range_data_vec_[i].second.returns) {
        Eigen::Matrix<T, 3, 1> world = transform * point.cast<T>();
        Eigen::Vector3f origin = Eigen::Vector3f(
            getScalar(transform.translation().x()),
            getScalar(transform.translation().y()),
            getScalar(transform.translation().z()));
        Eigen::Vector3f hit = Eigen::Vector3f(getScalar(world[0]),
                                              getScalar(world[1]),
                                              getScalar(world[2]));
        Eigen::Quaternionf orientation = Eigen::Quaternionf(
            getScalar(transform.rotation().w()),
            getScalar(transform.rotation().x()),
            getScalar(transform.rotation().y()),
            getScalar(transform.rotation().z()));

        std::vector<std::tuple<Eigen::Vector3f, double, double>> line =
            getCellInformationFromRayCast(origin, hit, counter);
        T prob_multiplicator = T(1.0);
        T sum = T(0.0);
        T cov_sum = T(0.0);
        double line_size_factor = 1.0 / line.size();
        double normalization_factor = 1.0 / (0.5 * hybrid_grid_->resolution())
            * std::exp(
                -1.0 / (0.5 * hybrid_grid_->resolution()) * 0.5
                    * hybrid_grid_->resolution());
        int incr = 0;
        Eigen::Vector3f point_in_world_base;
        for (std::tuple<Eigen::Vector3f, double, double>& tuple : line) {
          Eigen::Matrix<T, 3, 1> point_in_world;
          if (incr != line.size() - 1) {
            point_in_world_base = orientation.inverse()
                * (std::get<0>(tuple) * hybrid_grid_->resolution()) - origin;
            point_in_world = transform * point_in_world_base.cast<T>();
          } else {
            //point_in_world_base = orientation.inverse() * std::get<0>(tuple) - origin;
            point_in_world = world;

            //              LOG(INFO)<<"point:"<<point_in_world_base - point_cloud_[i];
            //              Eigen::Matrix<T, 3, 1> test = transform * point_in_world_base.cast<T>();
            //              Eigen::Vector3f test2 = Eigen::Vector3f(getScalar(test[0]),
            //                                                    getScalar(test[1]),
            //                                                    getScalar(test[2]));
            //              LOG(INFO)<<"point_after:"<<test2;
          }
          T prob;
          double factor;
          Eigen::Vector3f cell(getScalar(point_in_world[0]),
                               getScalar(point_in_world[1]),
                               getScalar(point_in_world[2]));
          Eigen::Array3i index = decay_grid_->GetCellIndex(cell);

//          if (decay_grid_->GetViewCount(index) > 0)
//          {
            T lambda = interpolated_decay_grid_.GetDecayRate(point_in_world[0],
                                                             point_in_world[1],
                                                             point_in_world[2]);
            //        Eigen::Array3i index(std::get<0>(tuple)[0], std::get<0>(tuple)[1], std::get<0>(tuple)[2]);
            //        T lambda = T(decay_grid_->GetDecayRate(index));
            prob = lambda * prob_multiplicator
                * ceres::exp(-lambda * std::get<1>(tuple)) * 1.0
                / normalization_factor;
            prob_multiplicator *= ceres::exp(-lambda * std::get<1>(tuple));

            //max value should be sqrt(3)/2 than the Factor should be 0
            //        if (prob > T(1.0))
            //          prob -= (prob - T(1.0));
//          }
//          else
//            prob = T(0.5);
          factor = -2 / sqrt(3) * std::get<2>(tuple) + 1;
          if (incr == line.size() - 1) {
            //Eigen::Matrix<T, 3, 1> diff = Eigen::Matrix<T, 3, 1> (world[0] - point_in_world[0], world[1] - point_in_world[1], world[2] - point_in_world[2]);
            //LOG(INFO)<<"x: "<<diff[0]<<"y: "<<diff[1]<<"z: "<<diff[3];
            //prob = interpolated_grid_.GetProbability(world[0], world[1], world[2]);
            //          LOG(INFO)<<getScalar(prob)<<" lambda "<<getScalar(lambda);
            sum += prob;
            sum = 1.0 - sum;
//            cov_sum += interpolated_decay_grid_.GetDecayRate(point_in_world[0], point_in_world[1], point_in_world[2])
//                            - T (getScalar(interpolated_decay_grid_.GetDecayRate(point_in_world[0], point_in_world[1], point_in_world[2])));
          } else {
            //          LOG(INFO)<<"free:"<<getScalar(prob)<<" lambda "<<getScalar(lambda);//<<" prob_mult:"<<prob_multiplicator;
            //          prob = interpolated_grid_.GetProbability(world[0], world[1], world[2]);
            sum += line_size_factor * prob;
          }
          incr++;
        }
//        T result = T(getScalar(sum)) + cov_sum;
        residual[counter] = scaling_factor_ * sum;
        counter++;

      }
    }
//      int counter = 0 ;
//      for (int i = 0; i < numSamplePoints; ++i) {
//    //    double point = minU
//    //        + (maxU - minU) / (double) (numSamplePoints - 1) * (double) i;
//        double point = ((common::ToSeconds(range_data_vec_[i].first - begin_)
//            / common::ToSeconds(end_ - begin_)));
//        std::array<T, output_dim> output;
//        nurbs.getPoint(&point, output.begin());
//        const transform::Rigid3<T> transform(
//            Eigen::Matrix<T, 3, 1>(output[0], output[1], output[2]),
//            Eigen::Quaternion<T>(output[6], output[3], output[4],
//                                 output[5]));
//        for (Eigen::Vector3f point : range_data_vec_[i].second.returns)
//        {
//          Eigen::Matrix<T, 3, 1> world = transform* point.cast<T>();
//          Eigen::Vector3f hit = Eigen::Vector3f(getScalar(world[0]),
//                                                getScalar(world[1]),
//                                                getScalar(world[2]));
//          Eigen::Array3i index = hybrid_grid_->GetCellIndex(hit);
//          T prob = T(hybrid_grid_->GetProbability(index));
//          const T probability =
//                      interpolated_grid_.GetProbability(world[0], world[1], world[2]);
//                  residual[counter] = scaling_factor_ * (1. - prob);
//          counter++;
//        }
//      }
//    int number_of_residuals = 0;
//    for (const std::pair<common::Time, sensor::RangeData>& data : range_data_vec_)
//      number_of_residuals += data.second.returns.size();
//    for (int i = 0; i < number_of_residuals; i++)
//      residual[i] = 0;
    return true;
  }

  std::vector<std::tuple<Eigen::Vector3f, double, double>> getCellInformationFromRayCast(
      Eigen::Vector3f& origin_f, Eigen::Vector3f& hit, int& counter) const {
    Eigen::Array3i origin = hybrid_grid_->GetCellIndex(origin_f);
    Eigen::Array3i current_key = hybrid_grid_->GetCellIndex(hit);
    std::vector<std::tuple<Eigen::Vector3f, double, double>> line;
    Eigen::Array3i hit_cell = hybrid_grid_->GetCellIndex(hit);
    //LOG(INFO)<<"hit count: "<<std::get<1>(*(hybrid_grid->mutable_value(hit_cell)));
    Eigen::Array3i direction = origin - hit_cell;
    Eigen::Vector3f slope(direction[0], direction[1], direction[2]);
    double distance = 0.5 * calculateRayLengthInVoxel(slope, origin, hit_cell)
        * decay_grid_->resolution();
    Eigen::Vector3f point = GetRayPoint(
        slope, Eigen::Vector3f(hit_cell[0], hit_cell[1], hit_cell[2]),
        origin_f);
    double distance_to_ray = DistanceToRay(
        slope, Eigen::Vector3f(hit_cell[0], hit_cell[1], hit_cell[2]),
        origin_f);
    line.push_back(std::make_tuple(hit, distance, distance_to_ray));

    //    LOG(INFO) << "hit: " << hit_cell << std::endl
    //              << " origin: " << origin << std::endl
    //              << " direction: " << direction;
    int step[3];
    double tMax[3];
    double tDelta[3];

    for (unsigned int i = 0; i < 3; ++i) {
      // compute step direction
      if (direction(i) > 0.0)
        step[i] = 1;
      else if (direction(i) < 0.0)
        step[i] = -1;
      else
        step[i] = 0;

      // compute tMax, tDelta
      if (step[i] != 0) {
        // corner point of voxel (in direction of ray)
        double voxelBorder = (step[i] * 0.5) / direction(i);
        tMax[i] = voxelBorder;
        tDelta[i] = 1.f / fabs(direction(i));

      } else {
        tMax[i] = std::numeric_limits<double>::max();
        tDelta[i] = std::numeric_limits<double>::max();
      }
    }

    // Incremental phase
    // ---------------------------------------------------------
    bool done = false;

    while (!done) {
      unsigned int dim;
      // find minimum tMax:
      if (tMax[0] < tMax[1]) {
        if (tMax[0] < tMax[2])
          dim = 0;
        else
          dim = 2;
      } else {
        if (tMax[1] < tMax[2])
          dim = 1;
        else
          dim = 2;
      }
      current_key[dim] += step[dim];
      tMax[dim] += tDelta[dim];
      double dist = calculateRayLengthInVoxel(slope, origin, current_key)
          * decay_grid_->resolution();
      point = GetRayPoint(
          slope,
          Eigen::Vector3f(current_key[0], current_key[1], current_key[2]),
          origin_f);
      distance_to_ray = DistanceToRay(
          slope,
          Eigen::Vector3f(current_key[0], current_key[1], current_key[2]),
          origin_f);      // * decay_grid_->resolution();
      line.insert(line.begin(), std::make_tuple(point, dist, distance_to_ray));
      if ((current_key(0) == origin(0) && current_key(1) == origin(1)
          && current_key(2) == origin(2)) || line.size() == line_size_) {
        done = true;
      } else {
      }
    }
    // end while
    //getRayInformation(line, slope, origin, counter);
    return line;
  }

  void getRayInformation(std::vector<std::pair<Eigen::Array3i, double>> & line,
                         Eigen::Vector3f& slope, Eigen::Array3i& origin,
                         int& counter) const {
//    double line_size_factor = 1.0 / line.size();
//    T prob_multiplicator = T(1);
//    Eigen::Vector3f origin_f(origin[0], origin[1], origin[2]);
//    T sum = T(0);
//    int incr = 0;
//    for (const std::pair<Eigen::Array3i,double>& cell : line)
//    {
//     double dist = cell.second;
//     T lambda, prob;
//     std::tuple<uint16, uint16, double> values = decay_grid_->value(cell.first);
//     Eigen::Vector3f point = GetRayPoint(slope, Eigen::Vector3f(cell.first[0],
//                                        cell.first[1],
//                                        cell.first[2]),
//                                              origin_f);
//     lambda = interpolated_decay_grid_.GetDecayRate(point_t[0], point_t[1], point_t[2]);
////     if (lambda != decay_grid_->GetDecayRate(cell.first)){
////       LOG(INFO)<<"lambda_not_inter:"<<decay_grid_->GetDecayRate(cell.first);
////       LOG(INFO)<<"lambda:"<<lambda;
////     }
//     prob = lambda * prob_multiplicator * ceres::exp(-lambda * dist);
////     if (prob > 1)
////     {
////       prob = 0.9;
////     }
////     else if (prob < 0.1){
////       prob = 0.1;
////     }
//       prob_multiplicator *= ceres::exp(-lambda * dist);
//     double factor = -2/sqrt(3) * DistanceToRay(slope, Eigen::Vector3f(cell.first[0],
//                                          cell.first[1],
//                                          cell.first[2]),
//                                                origin_f) + 1;
//
//     if (incr == line.size() - 1){
////       residual[counter] =  (1.0 - prob);
////       LOG(INFO)<<residual[counter];
//       counter++;
//     }
//     else
//       //sum +=line_size_factor * factor * T((0.1 - prob));
//     incr++;
//    }
  }

  double calculateRayLengthInVoxel(Eigen::Vector3f& slope,
                                   Eigen::Array3i& origin,
                                   Eigen::Array3i& cell) const {
    double t;
    Eigen::Vector3f cell_f(cell[0], cell[1], cell[2]);
    double scaling_factor = 1;
    if (origin[0] == cell[0] && origin[1] == cell[1] && origin[2] == cell[2])
      scaling_factor = 0.5;
    Eigen::Vector3f origin_f(origin[0], origin[1], origin[2]);
    std::vector<Eigen::Vector3f> intersections;
    std::vector<Eigen::Vector3f> normals = { Eigen::Vector3f::UnitX(),
        Eigen::Vector3f::UnitY(), Eigen::Vector3f::UnitZ() };
    for (Eigen::Vector3f& normal : normals) {
      for (double i = -0.5; i < 1.5; i++) {
        Eigen::Vector3f point_diff = cell_f + i * normal - origin_f;
        double inner_product_slope_normal = slope.dot(normal);
        if (inner_product_slope_normal == 0.0)
          continue;
        t = point_diff.dot(normal) / inner_product_slope_normal;
        Eigen::Vector3f intersection = t * slope + origin_f;
        if (checkInsideVoxel(cell_f, intersection)) {
          intersections.push_back(intersection);
        }
        if (intersections.size() == 2)
          return scaling_factor * ((intersections[1] - intersections[0]).norm());

      }
    }
    return 0;
  }
  bool checkInsideVoxel(Eigen::Vector3f& cell, Eigen::Vector3f& point) const {
    Eigen::Vector3f diff = cell - point;
    if (diff.norm() < sqrt(3) / 2.0)
      return true;
    else
      return false;
  }

  //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
  double DistanceToRay(Eigen::Vector3f slope, Eigen::Vector3f cell,
                       Eigen::Vector3f origin) const {
    Eigen::Vector3f diff = origin - cell;
    slope.normalize();
    return (diff - (diff.dot(slope) * slope)).norm();
  }

  //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
  Eigen::Vector3f GetRayPoint(Eigen::Vector3f slope, Eigen::Vector3f cell,
                              Eigen::Vector3f origin) const {
    Eigen::Vector3f diff = origin - cell;
    slope.normalize();
    return (origin - (diff.dot(slope) * slope));
  }

 private:
  const double scaling_factor_;
  const InterpolatedDecayGrid interpolated_decay_grid_;
  const InterpolatedGrid interpolated_grid_;
  const HybridGrid* hybrid_grid_;
  const HybridDecayGrid* decay_grid_;
  const common::Time begin_;
  const common::Time end_;
  const std::vector<std::pair<common::Time, sensor::RangeData>> range_data_vec_;
  std::shared_ptr<RayTracer> ray_tracer_;
  int line_size_;
};
}
}
}

#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_CONTROL_POINTS_FREE_SPACE_COST_FUNCTOR_H_ */
