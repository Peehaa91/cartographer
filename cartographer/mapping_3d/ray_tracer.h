/*
 * ray_tracer.h
 *
 *  Created on: Sep 11, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_RAY_TRACER_H_
#define CARTOGRAPHER_MAPPING_3D_RAY_TRACER_H_

#include "cartographer/mapping_3d/hybrid_grid.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/sensor/range_data.h"

namespace cartographer {
namespace mapping_3d {

class RayTracer {
 public:
  RayTracer();
  RayTracer(const int& line_size);
//  RayTracer(const RayTracer&) = delete;
//  RayTracer& operator=(const RayTracer&) = delete;
  std::vector<Eigen::Array3i> getLine(const Eigen::Vector3f& origin, const Eigen::Vector3f& hit,
                                      const HybridGrid* hybrid_grid);
 private:
  int max_line_size_;
};
}
}

#endif /* CARTOGRAPHER_MAPPING_3D_RAY_TRACER_H_ */
