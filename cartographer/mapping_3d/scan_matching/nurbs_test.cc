/*
 * nurbs_test.cc
 *
 *  Created on: Jul 17, 2017
 *      Author: schnattinger
 */

#include "gmock/gmock.h"
#include "cartographer/mapping_3d/scan_matching/nurbs.h"
#include <glog/logging.h>
#include <random>

namespace cartographer {
namespace mapping_3d {

TEST(NurbsTest, NurbsTest) {

  const int input_dim = 1;
  const int output_dim = 2;
  const int degree = 3;
  const KnotType knot_type = KnotType::UNIFORM;
  const WeightType weight_type = WeightType::RATIONAL;
  std::mt19937 rng;
  std::uniform_real_distribution<> dist(0.0, 10.0);
  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurb = Nurbs<
      double, input_dim, output_dim, knot_type, weight_type>();
  NurbsKnots<double, input_dim, output_dim, knot_type> knots = NurbsKnots<
      double, input_dim, output_dim, knot_type>();
  NurbsWeights<double, input_dim, output_dim, weight_type> weights =
      NurbsWeights<double, input_dim, output_dim, weight_type>();
  boost::multi_array<std::array<double, output_dim>, input_dim> points(
      boost::extents[10]);
  for (int i = 0; i < 10; ++i) {
    points[i][0] = (double) i;  //x component
    points[i][1] = dist(rng);  //y component
  }

  LOG(INFO)<<"shape: "<<*(points.shape());
  knots.create(degree, points.shape());
  std::array<std::vector<double>, input_dim> knot_vec = knots.getKnot();
  for (int i = 0; i < knot_vec[0].size(); i++)
    LOG(INFO)<<"knot: "<<knot_vec[0][i];
    //LOG(INFO)<<"knots: "<<knots.getKnot();
  weights.create(points.shape());

//  weights.create(weights_vec);
//  weights_vec = weights.getWeights();
//  for (int i = 0; i < weights_vec.size(); i++) {
//    LOG(INFO)<<"weight: "<<weights_vec[i];
//  }
  nurb.init(degree, knots, weights, points);
  double minU = nurb.getMinU(0);
  double maxU = nurb.getMaxU(0);
  //nurb.getPoint(1,0);
  //LOG(INFO)<<"valid: "<<nurb.validKnotVector();
  points = nurb.getPoints();
  for (int i = 0; i < 20; ++i) {
    double point = minU + (maxU - minU) / (double) (20 - 1) * (double) i;

    std::array<double, 1> output;
    nurb.getPoint(&point, output.begin());
    LOG(INFO)<<"x: "<<output[0]<<"  y: "<<output[1];

  }
  boost::multi_array<double, input_dim>& weights_vec = nurb.getWeights().getWeights();
  for (int i = 0; i < weights_vec.size(); i++) {
    LOG(INFO)<<"weight: "<<weights_vec[i];
    weights_vec[i] = static_cast<double>(i) / static_cast<double>(weights_vec.size());
    LOG(INFO)<<"weight: "<<weights_vec[i];
  }
//  weights.create(weights_vec);
 // nurb.init(degree, knots, weights, points);
  points = nurb.getPoints();
    for (int i = 0; i < 20; ++i) {
      double point = minU + (maxU - minU) / (double) (20 - 1) * (double) i;

      std::array<double, 1> output;
      nurb.getPoint(&point, output.begin());
      LOG(INFO)<<"x: "<<output[0]<<"  y: "<<output[1];

    }
}
}
}

