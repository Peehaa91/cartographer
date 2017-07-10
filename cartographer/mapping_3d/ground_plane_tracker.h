
#ifndef CARTOGRAPHER_MAPPING_GROUND_PLANE_TRACKER_H_
#define CARTOGRAPHER_MAPPING_GROUND_PLANE_TRACKER_H_

#include "Eigen/Geometry"
#include "cartographer/common/time.h"

namespace cartographer {
namespace mapping_3d {

// Keeps track of the orientation
class GroundPlaneTracker {
 public:
  GroundPlaneTracker();
  GroundPlaneTracker(Eigen::Quaterniond& quaternion, common::Time time);

  void Advance(common::Time time);

  void AddGroundPlaneObservation(const Eigen::Vector4d& coefficients);


  // Query the current orientation estimate.
  Eigen::Quaterniond orientation() const { return orientation_; }

 private:
  common::Time time_;
  common::Time last_coefficients_time_;
  Eigen::Quaterniond orientation_;
};
}
}
#endif
