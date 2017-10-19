/*
 * interpolated_decay_grid_test.cc
 *
 *  Created on: Sep 22, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_INTERPOLATED_DECAY_GRID_TEST_CC_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_INTERPOLATED_DECAY_GRID_TEST_CC_
#include "cartographer/mapping_3d/scan_matching/interpolated_decay_grid.h"

#include "Eigen/Core"
#include "cartographer/mapping_3d/hybrid_grid.h"
#include "gtest/gtest.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {
namespace {

class InterpolatedDecayGridTest : public ::testing::Test {
 protected:
  InterpolatedDecayGridTest()
      : hybrid_grid_(0.1f), interpolated_grid_(hybrid_grid_) {
    for (const auto& point :
         {Eigen::Vector3d(-3.f, 2.f, 0.f), Eigen::Vector3d(-4.f, 2.f, 0.f),
          Eigen::Vector3d(-5.f, 2.f, 0.f), Eigen::Vector3d(-6.f, 2.f, 0.f),
          Eigen::Vector3d(-6.f, 3.f, 1.f), Eigen::Vector3d(-6.f, 4.f, 2.f),
          Eigen::Vector3d(-7.f, 3.f, 1.f)}) {
      Eigen::Array3i array_point(point[0], point[1], point[2]);
//      LOG(INFO)<<"decay rate before: "<<hybrid_grid_.GetDecayRate(array_point);
      hybrid_grid_.increaseHitCount(array_point);
      double dist = 0.5;
      hybrid_grid_.increaseRayAccumulation(array_point, dist);
      LOG(INFO)<<"grid:"<<hybrid_grid_.GetDecayRate(array_point);
      LOG(INFO)<<"interpolated:"<<interpolated_grid_.GetDecayRate(point[0]*0.1 + 0.05,point[1]*0.1,point[2]*0.1);
//      LOG(INFO)<<"decay rate: "<<hybrid_grid_.GetDecayRate(array_point);
    }
  }

  float GetHybridGridProbability(float x, float y, float z) const {
    return hybrid_grid_.GetProbability(
        hybrid_grid_.GetCellIndex(Eigen::Vector3f(x, y, z)));
  }

  HybridDecayGrid hybrid_grid_;
  InterpolatedDecayGrid interpolated_grid_;
};

TEST_F(InterpolatedDecayGridTest, InterpolatesGridPoints) {
  for (int z = -1.; z < 3.; z += 1) {
    for (int y = 1.; y < 5.; y += 1) {
      for (int x = -8.; x < -2.; x += 1) {
//        Eigen::Array3i point(x,y,z);
//        EXPECT_NEAR(hybrid_grid_.GetDecayRate(point),
//                    interpolated_grid_.GetDecayRate(x, y, z), 1e-6);
//        LOG(INFO)<<"cell:"<<interpolated_grid_.GetDecayRate(x, y, z);
//        LOG(INFO)<<interpolated_grid_.GetDecayRate(x + 0.25, y + 0.25, z);
      }
    }
  }
}


}  // namespace
}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_INTERPOLATED_DECAY_GRID_TEST_CC_ */
