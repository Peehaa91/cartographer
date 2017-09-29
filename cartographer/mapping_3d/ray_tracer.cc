/*
 * ray_tracer.cpp
 *
 *  Created on: Sep 11, 2017
 *      Author: schnattinger
 */

#include "cartographer/mapping_3d/ray_tracer.h"

namespace cartographer {
namespace mapping_3d {
RayTracer::RayTracer() :
    max_line_size_(0){
}
RayTracer::RayTracer(const int& line_size) :
    max_line_size_(line_size){
  // TODO Auto-generated constructor stub

}


std::vector<Eigen::Array3i> RayTracer::getLine(const Eigen::Vector3f& origin, const Eigen::Vector3f& hit,
                                               const HybridGrid* hybrid_grid) {

  Eigen::Array3i current_key =
  hybrid_grid->GetCellIndex(hit);
  Eigen::Array3i origin_cell = hybrid_grid->GetCellIndex(origin);
  std::vector<Eigen::Array3i> line;
  Eigen::Array3i hit_cell = hybrid_grid->GetCellIndex(hit);
//  hybrid_grid->increaseHitCount(hit_cell);
  double distance = 0.5;
//  hybrid_grid->increaseRayAccumulation(hit_cell, distance);

  Eigen::Array3i direction = origin_cell - hit_cell;
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
  line.push_back(hit_cell);
  bool done = false;

  while (!done) {
    unsigned int dim;
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

    line.insert(line.begin(), current_key);
    if ((current_key(0) == origin_cell(0) && current_key(1) == origin_cell(1) &&
        current_key(2) == origin_cell(2)) || line.size() == max_line_size_) {
      double dist = 1;
//      hybrid_grid->increaseRayAccumulation(current_key, dist);
      done = true;
    }
    else {
      double dist = 1;
//      hybrid_grid->increaseRayAccumulation(current_key, dist);
        }
      }
  return line;
}
}
}

