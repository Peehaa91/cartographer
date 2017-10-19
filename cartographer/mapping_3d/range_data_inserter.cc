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
#include <chrono>

namespace cartographer {
namespace mapping_3d {

namespace {

void InsertMissesIntoGrid(const std::vector<uint16>& miss_table,
                          const Eigen::Vector3f& origin,
                          const sensor::PointCloud& returns,
                          HybridGrid* hybrid_grid,
                          const int num_free_space_voxels) {
  const Eigen::Array3i origin_cell = hybrid_grid->GetCellIndex(origin);
  int counter = 0;
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
      counter++;
      const Eigen::Array3i miss_cell =
          origin_cell + delta * position / num_samples;
      hybrid_grid->ApplyLookupTable(miss_cell, miss_table);
      //(INFO)<<"x: "<<miss_cell(0)<<" y: "<<miss_cell(1)<<" z: "<<miss_cell(2);
    }
  }
//  LOG(INFO)<<"updated cells: "<<counter;
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
          mapping::Odds(options_.miss_probability()))),
      ray_tracer_(RayTracer(options.num_free_space_voxels())){}

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
//  hybrid_grid->FinishUpdate();
}

void RangeDataInserter::Insert(const sensor::RangeData& range_data,
                               HybridDecayGrid* hybrid_decay_grid, HybridGrid* hybrid_grid) const {
//  LOG(INFO)<<"Ray Tracing started";
//  RayTracingInsert(range_data, hybrid_decay_grid, hybrid_grid);
//  hybrid_grid->updateWithDecayGrid(*hybrid_decay_grid);
//  LOG(INFO)<<"Ray Traycing finished";
}

void RangeDataInserter::RayTracingInsert(const sensor::RangeData& range_data,
                                         HybridDecayGrid* hybrid_decay_grid, HybridGrid* hybrid_grid) const {
  /// ----------  see OcTreeBase::computeRayKeys  -----------

  // Initialization phase
  // -------------------------------------------------------
  uint line_size = options_.num_free_space_voxels();
  Eigen::Array3i origin = hybrid_decay_grid->GetCellIndex(range_data.origin);
  std::vector<std::vector<Eigen::Array3i>> lines;
  for (const Eigen::Vector3f& hit : range_data.returns) {
    Eigen::Array3i current_key =
	    hybrid_decay_grid->GetCellIndex(hit);
    std::vector<Eigen::Array3i> line;
    Eigen::Array3i hit_cell = hybrid_decay_grid->GetCellIndex(hit);
    hybrid_decay_grid->increaseHitCount(hit_cell);
    hybrid_decay_grid->increaseViewCount(hit_cell);
    //double distance = 0.5;
    //LOG(INFO)<<"hit count: "<<std::get<1>(*(hybrid_grid->mutable_value(hit_cell)));
    Eigen::Array3i direction = origin - hit_cell;
    Eigen::Vector3f slope(direction[0], direction[1], direction[2]);
    double distance = 0.5*calculateRayLengthInVoxel(slope, origin, hit_cell) * hybrid_grid->resolution() ;
    hybrid_decay_grid->increaseRayAccumulation(hit_cell, distance);

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
    line.push_back(hit_cell);
    bool done = false;

    while (!done) {
      unsigned int dim;
      //LOG(INFO) << "0: " << tMax[0] << " 1: " << tMax[1] << " 2: " << tMax[2];
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
//      LOG(INFO) << "before "
//                << " cell_x: " << current_key(0) << std::endl
//                << " cell_y: " << current_key(1) << std::endl
//                << " cell_z: " << current_key(2);
      current_key[dim] += step[dim];
      tMax[dim] += tDelta[dim];

      // generate world coords from key
      //Eigen::Array3i pos = current_key.cast<int>();

      // LOG(INFO)<<" cell_x: "<<current_key(0)<<std::endl<<" cell_y:
      // "<<current_key(1)
      //    <<std::endl<<" cell_z: "<<current_key(2);
      line.insert(line.begin(), current_key);
      double dist = calculateRayLengthInVoxel(slope, origin, current_key) * hybrid_grid->resolution();
      hybrid_decay_grid->increaseViewCount(current_key);
//      LOG(INFO)<<current_key;
//      LOG(INFO)<<"dist: "<<dist;
      if ((current_key(0) == origin(0) && current_key(1) == origin(1) &&
          current_key(2) == origin(2)) || line.size() == line_size) {
        //double dist = calculateRayLengthInVoxel;
        hybrid_decay_grid->increaseRayAccumulation(current_key, dist);
        done = true;
      }
      else {
        //double dist = 1;
        hybrid_decay_grid->increaseRayAccumulation(current_key, dist);
        //LOG(INFO)<<"dist: "<<dist;
      }
    }
    // end while
    lines.push_back(line);
  }
  //LOG(INFO)<<"intersct_done";
  setProbabilities(lines, hybrid_decay_grid, hybrid_grid );
//  updateProbabilities(lines, hybrid_decay_grid, hybrid_grid );
//  hybrid_grid->updateWithDecayGrid(*hybrid_decay_grid);
  //LOG(INFO)<<"update prob done";
}

void RangeDataInserter::setProbabilities(const std::vector<std::vector<Eigen::Array3i>>& lines, HybridDecayGrid* hybrid_decay_grid, HybridGrid* hybrid_grid) const
{
  for (const std::vector<Eigen::Array3i>& line : lines){
    for (const Eigen::Array3i& index : line){
//     LOG(INFO)<<"index: "<<index<<" prob: "<<hybrid_decay_grid->GetProbability(index);
     hybrid_grid->SetProbability(index, hybrid_decay_grid->GetProbability(index));
    }
  }
}

double RangeDataInserter::DistanceToRay(Eigen::Vector3f slope, Eigen::Vector3f& cell, Eigen::Vector3f& origin) const
{
  Eigen::Vector3f diff = origin - cell;
  slope.normalize();
  Eigen::Vector3f line_point = origin - (diff.dot(slope)*slope);
  return (diff - (diff.dot(slope)*slope)).norm();

}

double RangeDataInserter::calculateRayLengthInVoxel(Eigen::Vector3f& slope,
                                                    Eigen::Array3i& origin,
                                                    Eigen::Array3i& cell) const
{
  double t;
  Eigen::Vector3f cell_f(cell[0], cell[1], cell[2]);
  double scaling_factor = 1;
  if (origin[0] == cell[0] && origin[1] == cell[1] && origin[2] == cell[2])
    scaling_factor = 0.5;
  Eigen::Vector3f origin_f(origin[0], origin[1], origin[2]);
  DistanceToRay(slope, cell_f, origin_f);
  std::vector<Eigen::Vector3f> intersections;
  std::vector<Eigen::Vector3f> normals = {Eigen::Vector3f::UnitX(),
      Eigen::Vector3f::UnitY(), Eigen::Vector3f::UnitZ()};
  for (Eigen::Vector3f& normal : normals)
  {
    for (double i = -0.5; i < 1.5; i++)
    {
      Eigen::Vector3f point_diff = cell_f + i*normal - origin_f;
      double inner_product_slope_normal = slope.dot(normal);
      if (inner_product_slope_normal == 0.0)
        continue;
      t = point_diff.dot(normal)/inner_product_slope_normal;
      Eigen::Vector3f intersection = t*slope + origin_f;
      if (checkInsideVoxel(cell_f ,intersection))
      {
        intersections.push_back(intersection);
      }
      if (intersections.size() == 2)
        return scaling_factor*((intersections[1] - intersections[0]).norm());

    }
  }
  return 0;
}
bool RangeDataInserter::checkInsideVoxel(Eigen::Vector3f& cell, Eigen::Vector3f& point) const
{
  Eigen::Vector3f diff = cell - point;
  if (diff.norm() <= sqrt(3)/2.0)
    return true;
  else
    return false;
}

float SlowValueToProbability(const uint16 value) {
  CHECK_GE(value, 0);
  CHECK_LE(value, 32767);
  if (value == mapping::kUnknownProbabilityValue) {
    // Unknown cells have kMinProbability.
    return mapping::kMinProbability;
  }
  const float kScale = (mapping::kMaxProbability - mapping::kMinProbability) / 32766.f;
  return value * kScale + (mapping::kMinProbability - kScale);
}

const std::vector<float>* PrecomputeValueToProbability() {
  std::vector<float>* result = new std::vector<float>;
  // Repeat two times, so that both values with and without the update marker
  // can be converted to a probability.
  for (int repeat = 0; repeat != 2; ++repeat) {
    for (int value = 0; value != 32768; ++value) {
      result->push_back(SlowValueToProbability(value));
    }
  }
  return result;
}
const std::vector<float>* const kValueToProbability =
    PrecomputeValueToProbability();
void RangeDataInserter::updateProbabilities(const std::vector<std::vector<Eigen::Array3i>>& lines, HybridDecayGrid* hybrid_decay_grid, HybridGrid* hybrid_grid) const
{
  std::chrono::high_resolution_clock::time_point begin_time =
  std::chrono::high_resolution_clock::now();
  //LOG(INFO)<<"update prob called";
  int count = 0;
  for (const std::vector<Eigen::Array3i>& line : lines)
  {
    double prob_multiplicator = 1;
    uint counter = 0;
    for ( std::vector<Eigen::Array3i>::const_iterator it = line.begin(); it != line.end(); it++)
    {
      double dist, prob;
      if (counter == line.size() - 1)
      {
        dist = 0.3;
//        prob = options_.hit_probability();
//        LOG(INFO)<<"hit";
//        hybrid_decay_grid->ApplyLookupTable(*it, hit_table_);
      }
      else{
//        LOG(INFO)<<"miss: ";
//        prob = options_.miss_probability();
//        hybrid_decay_grid->ApplyLookupTable(*it, miss_table_);
        dist = 1;
      }
      std::tuple<uint16, uint16, double, uint16>* values = hybrid_decay_grid->mutable_value(*it);
//      LOG(INFO)<<"val before: "<<std::get<0>(*values);
      double lambda = std::get<1>(*values)/std::get<2>(*values);
      prob = lambda * prob_multiplicator * exp(-lambda * dist);
      if (prob > 1)
      {
//        LOG(ERROR)<<"prob: "<<prob<<" lambda: "<<lambda<<" prob_mult: "<<prob_multiplicator<<" dist: "<<dist;
        prob = 0.9;
      }
      else if (prob < 0.1){
//        LOG(ERROR)<<"prob: "<<prob<<" lambda: "<<lambda<<" prob_mult: "<<prob_multiplicator<<" dist: "<<dist;
        prob = 0.1;
      }
      counter++;
      prob_multiplicator *= exp(-lambda * dist);
//      std::vector<uint16> table(mapping::ComputeLookupTableToApplyOdds(
//          mapping::Odds(prob)));
//      LOG(INFO)<<"pos: "<<*it;
 //     LOG(INFO)<<"update prob: "<<prob;
//      if (prob < 0.1)
//        prob = 0.1;
      hybrid_decay_grid->updateProbability(*it, prob, kValueToProbability);
//      LOG(INFO)<<"prob:"<<hybrid_decay_grid->GetProbability(*it);
//      hybrid_decay_grid->SetProbability(*it, prob);
//      LOG(INFO)<<"index: "<<*it<<std::endl<<"prob: "<<hybrid_decay_grid->GetProbability(*it)<<" value: "<<std::get<0>(*values);
      count++;
      hybrid_grid->SetProbability(*it, hybrid_decay_grid->GetProbability(*it));
//      hybrid_grid->ApplyLookupTable(*it, table);
//      if (counter == line.size())
//      {
//        LOG(INFO)<<"prob_after: "<<hybrid_grid->GetProbability(index)<<" index: "<<index;
//      }

    }
//    LOG(INFO)<<"line number: "<<count<<" from:"<<lines.size()<<" finished";
//    count++;
//    hybrid_grid->FinishUpdate();
//    LOG(INFO)<<"finished called";
  }
  auto duration = std::chrono::high_resolution_clock::now() - begin_time;
  auto delta_t = std::chrono::duration_cast<std::chrono::microseconds>(
      duration);
  LOG(INFO)<<"prob update time: "<<delta_t.count();
  LOG(INFO)<<"counter: "<<count;
}
}  // namespace mapping_3d
}  // namespace cartographer
