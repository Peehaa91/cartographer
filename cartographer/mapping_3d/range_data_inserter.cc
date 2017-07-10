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

#include "cartographer/mapping_3d/range_data_inserter.h"

#include "Eigen/Core"
#include "cartographer/mapping/probability_values.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping_3d {

namespace {

void InsertMissesIntoGrid(const std::vector<uint16>& miss_table,
                          const Eigen::Vector3f& origin,
                          const sensor::PointCloud& returns,
                          HybridGrid* hybrid_grid,
                          const int num_free_space_voxels) {
  const Eigen::Array3i origin_cell = hybrid_grid->GetCellIndex(origin);
  for (const Eigen::Vector3f& hit : returns) {
    const Eigen::Array3i hit_cell = hybrid_grid->GetCellIndex(hit);

    const Eigen::Array3i delta = hit_cell - origin_cell;
    const int num_samples = delta.cwiseAbs().maxCoeff();
    CHECK_LT(num_samples, 1 << 15);
    // 'num_samples' is the number of samples we equi-distantly place on the
    // line between 'origin' and 'hit'. (including a fractional part for sub-
    // voxels) It is chosen so that between two samples we change from one voxel
    // to the next on the fastest changing dimension.
    //
    // Only the last 'num_free_space_voxels' are updated for performance.
    for (int position = std::max(0, num_samples - num_free_space_voxels);
         position < num_samples; ++position) {
      const Eigen::Array3i miss_cell =
          origin_cell + delta * position / num_samples;
      hybrid_grid->ApplyLookupTable(miss_cell, miss_table);
    }
  }
}

}  // namespace

proto::RangeDataInserterOptions CreateRangeDataInserterOptions(
    common::LuaParameterDictionary* parameter_dictionary) {
  proto::RangeDataInserterOptions options;
  options.set_hit_probability(
      parameter_dictionary->GetDouble("hit_probability"));
  options.set_miss_probability(
      parameter_dictionary->GetDouble("miss_probability"));
  options.set_num_free_space_voxels(
      parameter_dictionary->GetInt("num_free_space_voxels"));
  CHECK_GT(options.hit_probability(), 0.5);
  CHECK_LT(options.miss_probability(), 0.5);
  return options;
}

RangeDataInserter::RangeDataInserter(
    const proto::RangeDataInserterOptions& options)
    : options_(options),
      hit_table_(mapping::ComputeLookupTableToApplyOdds(
          mapping::Odds(options_.hit_probability()))),
      miss_table_(mapping::ComputeLookupTableToApplyOdds(
          mapping::Odds(options_.miss_probability()))) {}

void RangeDataInserter::Insert(const sensor::RangeData& range_data,
                               HybridGrid* hybrid_grid) const {
  CHECK_NOTNULL(hybrid_grid)->StartUpdate();

  for (const Eigen::Vector3f& hit : range_data.returns) {
    const Eigen::Array3i hit_cell = hybrid_grid->GetCellIndex(hit);
    hybrid_grid->ApplyLookupTable(hit_cell, hit_table_);
  }

  // By not starting a new update after hits are inserted, we give hits priority
  // (i.e. no hits will be ignored because of a miss in the same cell).
  InsertMissesIntoGrid(miss_table_, range_data.origin, range_data.returns,
                       hybrid_grid, options_.num_free_space_voxels());
}

void RangeDataInserter::Insert(const sensor::RangeData& range_data,
                               HybridDecayGrid* hybrid_grid) const {
  RayTracingInsert(range_data, hybrid_grid);
}

void RangeDataInserter::RayTracingInsert(const sensor::RangeData& range_data,
                                                         HybridDecayGrid* hybrid_grid) const{

  /// ----------  see OcTreeBase::computeRayKeys  -----------

    // Initialization phase -------------------------------------------------------
    Eigen::Array3f current_key = hybrid_grid->GetCellIndex(range_data.origin).cast<float>();
    std::vector<std::vector<Eigen::Array3i>> lines;
    for (const Eigen::Vector3f& hit : range_data.returns)
    {
      std::vector<Eigen::Array3i> line;
      Eigen::Array3i hit_cell = hybrid_grid->GetCellIndex(hit);
      Eigen::Vector3f direction = hit - range_data.origin;
      LOG(INFO)<<"hit: "<<hit_cell<<std::endl<<" origin: "<<current_key<<std::endl<<" direction: "<<direction;
      int step[3];
      double tMax[3];
      double tDelta[3];

      for(unsigned int i=0; i < 3; ++i) {
        // compute step direction
        if (direction(i) > 0.0) step[i] =  1;
        else if (direction(i) < 0.0)   step[i] = -1;
        else step[i] = 0;

        // compute tMax, tDelta
        if (step[i] != 0) {
          // corner point of voxel (in direction of ray)
          double voxelBorder = current_key[i];
          voxelBorder += double(step[i] * hybrid_grid->resolution() * 0.5);

          tMax[i] = ( voxelBorder - range_data.origin(i) ) / direction(i);
          tDelta[i] = hybrid_grid->resolution() / fabs( direction(i) );
        }
        else {
          tMax[i] =  std::numeric_limits<double>::max();
          tDelta[i] = std::numeric_limits<double>::max();
        }
      }


      // Incremental phase  ---------------------------------------------------------

      bool done = false;

      while (!done) {
        unsigned int dim;
        LOG(INFO)<<"0: "<<tMax[0]<<" 1: "<<tMax[1]<<" 2: "<<tMax[2];
        // find minimum tMax:
        if (tMax[0] < tMax[1]){
          if (tMax[0] < tMax[2]) dim = 0;
          else                   dim = 2;
        }
        else {
          if (tMax[1] < tMax[2]) dim = 1;
          else                   dim = 2;
        }

        // advance in direction "dim"
        LOG(INFO)<<"before "<<" cell_x: "<<current_key(0)<<std::endl<<" cell_y: "<<current_key(1)
            <<std::endl<<" cell_z: "<<current_key(2);
        current_key[dim] += step[dim];
        tMax[dim] += tDelta[dim];


        // generate world coords from key
        Eigen::Array3i pos = current_key.cast<int>();
        LOG(INFO)<<" cell_x: "<<current_key(0)<<std::endl<<" cell_y: "<<current_key(1)
            <<std::endl<<" cell_z: "<<current_key(2);
        line.push_back(pos);
        if (pos(0) == hit_cell(0) && pos(1) == hit_cell(1) && pos(2) == hit_cell(2))
          done = true;
      }
      // end while
      lines.push_back(line);
    }
}

}  // namespace mapping_3d
}  // namespace cartographer
