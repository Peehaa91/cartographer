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
    // LOG(INFO)<<"origin: "<<origin_cell;
    for (int position = std::max(0, num_samples - num_free_space_voxels);
         position < num_samples; ++position) {
      const Eigen::Array3i miss_cell =
          origin_cell + delta * position / num_samples;
      hybrid_grid->ApplyLookupTable(miss_cell, miss_table);
      //(INFO)<<"x: "<<miss_cell(0)<<" y: "<<miss_cell(1)<<" z: "<<miss_cell(2);
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
  CHECK_NOTNULL(hybrid_grid);

  for (const Eigen::Vector3f& hit : range_data.returns) {
    const Eigen::Array3i hit_cell = hybrid_grid->GetCellIndex(hit);
    hybrid_grid->ApplyLookupTable(hit_cell, hit_table_);
  }

  // By not starting a new update after hits are inserted, we give hits priority
  // (i.e. no hits will be ignored because of a miss in the same cell).
  InsertMissesIntoGrid(miss_table_, range_data.origin, range_data.returns,
                       hybrid_grid, options_.num_free_space_voxels());
  hybrid_grid->FinishUpdate();
}

void RangeDataInserter::Insert(const sensor::RangeData& range_data,
                               HybridDecayGrid* hybrid_grid) const {
  RayTracingInsert(range_data, hybrid_grid);
}

void RangeDataInserter::RayTracingInsert(const sensor::RangeData& range_data,
                                         HybridDecayGrid* hybrid_grid) const {
  /// ----------  see OcTreeBase::computeRayKeys  -----------

  // Initialization phase
  // -------------------------------------------------------
  Eigen::Array3i origin = hybrid_grid->GetCellIndex(range_data.origin);
  std::vector<std::vector<Eigen::Array3i>> lines;
  for (const Eigen::Vector3f& hit : range_data.returns) {
    Eigen::Array3f current_key =
	    hybrid_grid->GetCellIndex(range_data.origin).cast<float>();
    std::vector<Eigen::Array3i> line;
    Eigen::Array3i hit_cell = hybrid_grid->GetCellIndex(hit);
    hybrid_grid->increaseHitCount(hit_cell);
    LOG(INFO)<<"hit count: "<<std::get<1>(*(hybrid_grid->mutable_value(hit_cell)));
    Eigen::Array3i direction = hit_cell - origin;
    LOG(INFO) << "hit: " << hit_cell << std::endl
              << " origin: " << origin << std::endl
              << " direction: " << direction;
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
      LOG(INFO) << "0: " << tMax[0] << " 1: " << tMax[1] << " 2: " << tMax[2];
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

      // advance in direction "dim"
      LOG(INFO) << "before "
                << " cell_x: " << current_key(0) << std::endl
                << " cell_y: " << current_key(1) << std::endl
                << " cell_z: " << current_key(2);
      current_key[dim] += step[dim];
      tMax[dim] += tDelta[dim];

      // generate world coords from key
      Eigen::Array3i pos = current_key.cast<int>();

      // LOG(INFO)<<" cell_x: "<<current_key(0)<<std::endl<<" cell_y:
      // "<<current_key(1)
      //    <<std::endl<<" cell_z: "<<current_key(2);
      line.push_back(pos);
      if (pos(0) == hit_cell(0) && pos(1) == hit_cell(1) &&
          pos(2) == hit_cell(2)) {
        double dist = 0.5;
        hybrid_grid->increaseRayAccumulation(pos, dist);
        dist = std::get<2>(*(hybrid_grid->mutable_value(pos)));
    	  LOG(INFO)<<"done";
        done = true;
      }
      else {
        double dist = 1;
        hybrid_grid->increaseRayAccumulation(pos, dist);
        dist = std::get<2>(*(hybrid_grid->mutable_value(pos)));
        LOG(INFO)<<"dist: "<<dist;
      }
    }
    // end while
    lines.push_back(line);
  }
  updateProbabilities(lines, hybrid_grid);
}
void RangeDataInserter::updateProbabilities(const std::vector<std::vector<Eigen::Array3i>>& lines, HybridDecayGrid* hybrid_grid) const
{
  for (const std::vector<Eigen::Array3i>& line : lines)
  {
    double prob_multiplicator = 1;
    int counter = 0;
    for (const Eigen::Array3i& index : line)
    {
      double dist;
      if (counter == line.size() - 1)
      {
        dist = 0.5;
        LOG(INFO)<<"prob_before: "<<hybrid_grid->GetProbability(index);
      }
      else
        dist = 1;
      std::tuple<uint16, uint16, double>* values = hybrid_grid->mutable_value(index);
      double lambda = std::get<1>(*values)/std::get<2>(*values);
      double prob = lambda * prob_multiplicator * exp(-lambda * dist);
      counter++;
      prob_multiplicator *= exp(-lambda * dist);
      std::vector<uint16> table(mapping::ComputeLookupTableToApplyOdds(
          mapping::Odds(prob)));
      hybrid_grid->ApplyLookupTable(index, table);
      if (counter == line.size())
      {
        LOG(INFO)<<"prob_after: "<<hybrid_grid->GetProbability(index)<<" index: "<<index;
      }

    }
    hybrid_grid->FinishUpdate();
  }
}
}  // namespace mapping_3d
}  // namespace cartographer
