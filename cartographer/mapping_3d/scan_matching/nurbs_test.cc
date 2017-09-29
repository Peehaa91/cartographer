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
  const WeightType weight_type = WeightType::NON_RATIONAL;
  std::mt19937 rng;
  std::uniform_real_distribution<> dist(0.0, 10.0);
  Nurbs<double, input_dim, output_dim, knot_type, weight_type> nurb = Nurbs<
      double, input_dim, output_dim, knot_type, weight_type>();
  NurbsKnots<double, input_dim, output_dim, knot_type> knots = NurbsKnots<
      double, input_dim, output_dim, knot_type>();
  NurbsWeights<double, input_dim, output_dim, weight_type> weights =
      NurbsWeights<double, input_dim, output_dim, weight_type>();
  boost::multi_array<std::array<double, output_dim>, input_dim> points(
      boost::extents[5]);
  for (int i = 0; i < 5; ++i) {
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
//  boost::multi_array<double, input_dim>& weights_vec = nurb.getWeights().getWeights();
//  for (int i = 0; i < weights_vec.size(); i++) {
//    LOG(INFO)<<"weight: "<<weights_vec[i];
//    weights_vec[i] = 0.1;
//    LOG(INFO)<<"weight: "<<weights_vec[i];
//  }
//  weights.create(weights_vec);
 // nurb.init(degree, knots, weights, points);
  points = nurb.getPoints();
  for (int i = 0; i < 20; ++i) {
    double point = minU + (maxU - minU) / (double) (20 - 1) * (double) i;

    std::array<double, 1> output;
    nurb.getPoint(&point, output.begin());
    LOG(INFO)<<"x: "<<output[0]<<"  y: "<<output[1];

  }
  Nurbs<double, input_dim, output_dim, KnotType::NON_UNIFORM, weight_type> nurb_deriv = Nurbs<
      double, input_dim, output_dim, KnotType::NON_UNIFORM, weight_type>();
  NurbsKnots<double, input_dim, output_dim, KnotType::NON_UNIFORM> knots_deriv = NurbsKnots<
      double, input_dim, output_dim, KnotType::NON_UNIFORM>();
  NurbsWeights<double, input_dim, output_dim, weight_type> weights_deriv =
      NurbsWeights<double, input_dim, output_dim, weight_type>();
  boost::multi_array<std::array<double, output_dim>, input_dim> points_deriv(
      boost::extents[4]);
  for (int i = 0; i < 4; ++i) {
    points_deriv[i][0] = degree/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][0] - points[i][0]);
    points_deriv[i][1] = degree/(knot_vec[0][i + degree + 1] - knot_vec[0][i + 1]) * (points[i + 1][1] - points[i][1]);

  }
  knots_deriv.create(degree - 1, points_deriv.shape());
  std::array<std::vector<double>, input_dim> knot_vec_2 = std::array<std::vector<double>, input_dim>();
  const std::array<std::vector<double>, input_dim>& knot_vec_init = knots_deriv.getKnot();
  LOG(INFO)<<"size old knots"<<knot_vec_init[0].size();
  for (int i = 0; i < knot_vec_init[0].size(); i++){
    LOG(INFO)<<"knots before: "<<knot_vec_init[0][i];
    knot_vec_2[0].push_back(knot_vec[0][i+1]);
  }
  LOG(INFO)<<"size new knots"<<knot_vec_2[0].size();
  knots_deriv.setKnots(knot_vec_2);
  const std::array<std::vector<double>, input_dim>& knot_vec_3 = knots_deriv.getKnot();
  for (int i = 0; i < knot_vec_3[0].size(); i++){
    LOG(INFO)<<"knot: "<<knot_vec_3[0][i];
    //knot_vec_2[0][i] = knot_vec[0][i+1];
  }

  weights_deriv.create(points_deriv.shape());
  nurb_deriv.init(degree -1, knots_deriv, weights_deriv, points_deriv);
  minU = nurb_deriv.getMinU(0);
  maxU = nurb_deriv.getMaxU(0);
  //nurb.getPoint(1,0);
  //LOG(INFO)<<"valid: "<<nurb.validKnotVector();
  points = nurb_deriv.getPoints();
  for (int i = 0; i < 20; ++i) {
    double point = minU + (maxU - minU) / (double) (20 - 1) * (double) i;

    std::array<double, 1> output;
    nurb_deriv.getPoint(&point, output.begin());
    LOG(INFO)<<"x: "<<output[0]<<"  y: "<<output[1];
  }



}
}
}

