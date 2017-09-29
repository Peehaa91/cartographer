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

#ifndef CARTOGRAPHER_MAPPING_3D_RANGE_DATA_INSERTER_H_
#define CARTOGRAPHER_MAPPING_3D_RANGE_DATA_INSERTER_H_

#include <math.h>

#include "cartographer/mapping_3d/hybrid_grid.h"
#include "cartographer/mapping_3d/proto/range_data_inserter_options.pb.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/mapping_3d/ray_tracer.h"

namespace cartographer {
namespace mapping_3d {

proto::RangeDataInserterOptions CreateRangeDataInserterOptions(
    common::LuaParameterDictionary* parameter_dictionary);

class RangeDataInserter {
 public:
  explicit RangeDataInserter(const proto::RangeDataInserterOptions& options);

  RangeDataInserter(const RangeDataInserter&) = delete;
  RangeDataInserter& operator=(const RangeDataInserter&) = delete;

  // Inserts 'range_data' into 'hybrid_grid'.
  void Insert(const sensor::RangeData& range_data,
              HybridGrid* hybrid_grid) const;

  // Inserts 'range_data' into 'hybrid_grid'.
  void Insert(const sensor::RangeData& range_data,
              HybridDecayGrid* hybrid_decay_grid, HybridGrid* hybrid_grid) const;

  // Inserts 'range data' via ray tracing into hybrid_grid
  void RayTracingInsert(const sensor::RangeData& range_data,
		  	  HybridDecayGrid* hybrid_decay_grid, HybridGrid* hybrid_grid) const;

 private:
  void updateProbabilities(const std::vector<std::vector<Eigen::Array3i>>& lines,
                           HybridDecayGrid* hybrid_decay_grid,
                           HybridGrid* hybrid_grid) const;

  void setProbabilities(const std::vector<std::vector<Eigen::Array3i>>& lines,
                        HybridDecayGrid* hybrid_decay_grid,
                        HybridGrid* hybrid_grid) const;

  double calculateRayLengthInVoxel(Eigen::Vector3f& slope,
                                   Eigen::Array3i& origin,
                                   Eigen::Array3i& cell) const;

  bool checkInsideVoxel(Eigen::Vector3f& cell, Eigen::Vector3f& point) const;

  double DistanceToRay(Eigen::Vector3f slope, Eigen::Vector3f& cell, Eigen::Vector3f& origin) const;
  const proto::RangeDataInserterOptions options_;
  const std::vector<uint16> hit_table_;
  const std::vector<uint16> miss_table_;
  RayTracer ray_tracer_;
};

}  // namespace mapping_3d
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_3D_RANGE_DATA_INSERTER_H_
