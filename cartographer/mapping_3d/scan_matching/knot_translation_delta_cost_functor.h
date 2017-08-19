/*
 * knot_translation_delta_cost_functor.h
 *
 *  Created on: Aug 18, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_KNOT_TRANSLATION_DELTA_COST_FUNCTOR_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_KNOT_TRANSLATION_DELTA_COST_FUNCTOR_H_


#include "Eigen/Core"
#include "cartographer/transform/rigid_transform.h"

namespace cartographer {
namespace mapping_3d {
namespace scan_matching {

// Computes the cost of translating the initial pose estimate. Cost increases
// with the solution's distance from the initial estimate.
class KnotTranslationDeltaCostFunctor {
 public:
  // Constructs a new TranslationDeltaCostFunctor from the given
  // 'initial_pose_estimate'.
  explicit KnotTranslationDeltaCostFunctor(
      const double scaling_factor,
      std::array<double, 6> pose)
      : scaling_factor_(scaling_factor),
        x_(pose[0]),
        y_(pose[1]),
        z_(pose[2]) {}

  KnotTranslationDeltaCostFunctor(const KnotTranslationDeltaCostFunctor&) = delete;
  KnotTranslationDeltaCostFunctor& operator=(const KnotTranslationDeltaCostFunctor&) =
      delete;

  template <typename T>
  bool operator()(const T* const translation, T* residual) const {
    residual[0] = scaling_factor_ * (translation[0] - x_);
    residual[1] = scaling_factor_ * (translation[1] - y_);
    residual[2] = scaling_factor_ * (translation[2] - z_);
    return true;
  }

 private:
  const double scaling_factor_;
  const double x_;
  const double y_;
  const double z_;
};

}  // namespace scan_matching
}  // namespace mapping_3d
}  // namespace cartographer


#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_KNOT_TRANSLATION_DELTA_COST_FUNCTOR_H_ */
