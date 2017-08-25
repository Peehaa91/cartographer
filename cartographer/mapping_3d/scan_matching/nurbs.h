/*
 * nurbs.h
 *
 *  Created on: Jul 16, 2017
 *      Author: schnattinger
 */

#ifndef CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_H_
#define CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_H_
#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <array>
#include <numeric>
#include <limits>

#include <boost/multi_array.hpp>

namespace cartographer{
namespace mapping_3d{
enum class KnotType {
  UNIFORM,
  NON_UNIFORM
};

enum class WeightType {
  NON_RATIONAL,
  RATIONAL
};

//forward declaration
namespace NurbsCeresPrivate {
  template <std::size_t inputDim, WeightType weightType>
  struct CreateParameterVector;
}

  template <typename Iter>
  bool isBaseCorrect(Iter begin, Iter end) {
    typedef typename std::iterator_traits<Iter>::value_type ValueT;

    bool error = false;
    std::for_each(begin, end, [&error] (ValueT val) {
      if (val < ValueT(-1e7)) {
        error = true;
      }
    });
    if (error == true) {
      return false;
    }

    auto val = std::accumulate(begin, end, ValueT(0.0));

    using std::abs;

    auto epsilon = std::numeric_limits<ValueT>::epsilon();
    return (abs(val - ValueT(1.0)) < ValueT(10 *  std::distance(begin, end)) * epsilon);
  }

  template <typename T, size_t NumDims, typename Allocator>
  std::array<int, NumDims> getSizeVector(const boost::multi_array<T, NumDims, Allocator>& multiArray) {
    auto size = multiArray.shape();
    std::array<int, NumDims> sizes;
    std::transform(size, size + NumDims, sizes.begin(), [] (size_t val) {
      return (int)val;
    });
    return sizes;
  }

  template <typename T, size_t NumDims, typename Allocator1, typename Allocator2>
  void assignMultiArray(const boost::multi_array<T, NumDims, Allocator1>& from, boost::multi_array<T, NumDims, Allocator2>& to) {
    to.resize(getSizeVector(from));
    to = from;
  }

  template <typename KnotIterT>
  void createUniformKnot(int numPoint, int degree, KnotIterT knotIter) {
    assert(numPoint >= degree - 1);
    typedef typename std::remove_reference<decltype(*knotIter)>::type value_type;
    std::fill(knotIter, knotIter + degree, (value_type)0.0);

    knotIter = knotIter + degree;
    auto totalSteps = numPoint - degree + 1;
    for (int i = 0; i < totalSteps; ++i, ++knotIter) {
      *knotIter = (value_type)i / (value_type)(totalSteps - 1);
    }

    std::fill(knotIter, knotIter + degree, (value_type)1.0);
  }

  template <typename T, std::size_t inputDim, std::size_t curDim>
  struct getPoint {
    template <typename Func>
    static void apply(const std::array<std::vector<T>, inputDim>& ns, const std::array<int, inputDim>& ks,
              const std::array<int, inputDim>& degrees, std::array<int, inputDim>& weightIdx, T nMul, Func&& func) {
      const auto& p = degrees[curDim];
      const auto& k = ks[curDim];
      const auto& n = ns[curDim];
      auto& curWeightIdx = weightIdx[curDim];

      for (int i = 0; i <= p; ++i) {
        curWeightIdx = k - i;
        getPoint<T, inputDim, curDim + 1>::apply(ns, ks, degrees, weightIdx, nMul * n[i], std::forward<Func>(func));
      }
    }
  };

  template <typename T, std::size_t inputDim>
  struct getPoint < T, inputDim, inputDim > {
    template <typename Func>
    static void apply(const std::array<std::vector<T>, inputDim>& /*ns*/, const std::array<int, inputDim>& /*ks*/,
              const std::array<int, inputDim>& /*degrees*/, std::array<int, inputDim>& weightIdx, T nMul, Func&& func) {

      func(weightIdx, nMul);
    }
  };

  template <typename T, std::size_t inputDim, typename Func>
  void applyGetPoint(const std::array<std::vector<T>, inputDim>& ns, const std::array<int, inputDim>& ks,
             const std::array<int, inputDim>& degrees, Func&& func) {

    std::array<int, inputDim> weightIdx; //not needed to be initialized
    getPoint<T, inputDim, 0>::apply(ns, ks, degrees, weightIdx, T(1.0), std::forward<Func>(func));
  }


template <typename T, std::size_t inputDim, std::size_t outputDim, KnotType knotType>
class NurbsKnots {
};

template <typename T, std::size_t inputDim, std::size_t outputDim>
class NurbsKnots < T, inputDim, outputDim, KnotType::NON_UNIFORM > {
private:
  typedef std::array<std::vector<T>, inputDim> KnotsT;

public:
  NurbsKnots() = default;

  NurbsKnots(const NurbsKnots<T, inputDim, outputDim, KnotType::NON_UNIFORM>& o) : knots_(o.getKnot()) {
  }
  NurbsKnots(const NurbsKnots<T, inputDim, outputDim, KnotType::UNIFORM>& o) : knots_(o.getKnot()) {
  }

  NurbsKnots<T, inputDim, outputDim, KnotType::NON_UNIFORM>& operator=(const NurbsKnots<T, inputDim, outputDim, KnotType::NON_UNIFORM>& o) {
    if (&o != this) {
      knots_ = o.getKnot();
    }

    return *this;
  }
  NurbsKnots<T, inputDim, outputDim, KnotType::NON_UNIFORM>& operator=(const NurbsKnots<T, inputDim, outputDim, KnotType::UNIFORM>& o) {
    knots_ = o.getKnot();
    return *this;
  }

  void create(const KnotsT& knots) {
    knots_ = knots;
  }

  template <typename DegreesIter, typename PointsIter>
  void create(DegreesIter degreesIter, PointsIter numPointsIter) {
    for (int i = 0; i < (int)inputDim; ++i) {
      const auto& numPoint = *numPointsIter;

      const auto& degree = *degreesIter;
      auto& knot = knots_[i];

      knot.resize(numPoint + degree + 1);
      createUniformKnot(numPoint, degree, knot.begin());

      ++numPointsIter;
      ++degreesIter;
    }
  }

  template <typename PointsIter>
  void create(int degree, PointsIter numPointsIter) {
    for (int i = 0; i < (int)inputDim; ++i) {
      const auto& numPoint = *numPointsIter;

      auto& knot = knots_[i];

      knot.resize(numPoint + degree + 1);
      createUniformKnot((int)numPoint, degree, knot.begin());

      ++numPointsIter;
    }
  }

  const KnotsT& getKnot() const {
    return knots_;
  }

  int getKnotIdx(int dim, int degree, T u) const {
    const auto& knot = knots_[dim];
    return (int)std::distance(knot.begin(),
                  std::lower_bound(knot.begin() + degree + 1, knot.end() - degree - 1, u)) - 1;
  }

  template <typename OutIter>
  void computeBasis(int k, int dim, int degree, T u, OutIter n) const {
    const auto& knot = knots_[dim];
    n[0] = T(1.0);

    for (int i = 1; i <= degree; ++i) {
      n[i] = T(0.0);

      for (int j = i - 1; j >= 0; --j) {
        T a = T(0.0);
        if (knot[k + i - j] != knot[k - j]) {
          a = (u - knot[k - j]) / (knot[k + i - j] - knot[k - j]);
        }

        n[j + 1] += n[j] * (T(1.0) - a);
        n[j] = n[j] * a;
      }
    }
  }

  T getMinU(int inputD) const {
    //TODO: implement
    //return knots_[inputD][degrees_[inputD]];
    return T(0.0);
  }

  T getMaxU(int inputD) const {
    //TODO: implement
    //return knots_[inputD][knots_[inputD].size() - degrees_[inputD] - 1];
    return T(1.0);
  }

  template <typename Archive>
  void serialize(Archive& archive) {
    archive(knots_);
  }

private:
  KnotsT knots_; //one knot vector for each input dimension
};

template <typename T, std::size_t inputDim, std::size_t outputDim>
class NurbsKnots < T, inputDim, outputDim, KnotType::UNIFORM > {
private:
  typedef std::array<std::vector<T>, inputDim> KnotsT;

public:
  KnotsT getKnot() const {
    std::array<std::vector<T>, inputDim> knots;
    for (int dim = 0; dim < (int)inputDim; ++dim) {
      knots[dim].resize(pointDims_[dim] + degrees_[dim] + 1);
      createUniformKnot(pointDims_[dim], degrees_[dim], knots[dim].begin());
    }

    return knots;
  }

  template <typename PointsIter>
  void create(int degree, PointsIter numPointsBeginIter) {
    std::array<int, inputDim> degrees;
    std::fill(degrees.begin(), degrees.end(), degree);
    create(degrees.begin(), numPointsBeginIter);
  }

  template <typename DegreesIter, typename PointsIter>
  void create(DegreesIter degreesIter, PointsIter numPointsBeginIter) {
    auto numPointsEndIter = numPointsBeginIter;
    std::advance(numPointsEndIter, inputDim);

    std::copy_n(degreesIter, inputDim, degrees_.begin());

    std::transform(numPointsBeginIter, numPointsEndIter, pointDims_.begin(), [] (const typename std::iterator_traits<PointsIter>::reference val) -> int {
      return static_cast<int>(val);
    });

    std::transform(numPointsBeginIter, numPointsEndIter, degreesIter, knotDeltas_.begin(), [] (
      const typename std::iterator_traits<PointsIter>::reference pointDim,
      const typename std::iterator_traits<DegreesIter>::reference degree) {

      return (1.0) / (pointDim - degree);
    });
    for (int i = 0; i < (int)inputDim; ++i) {
      const auto& numPoint = *numPointsBeginIter;

      const auto& degree = *degreesIter;
      auto& knot = knots_[i];

      knot.resize(numPoint + degree + 1);
      createUniformKnot(static_cast<int>(numPoint), degree, knot.begin());

      ++numPointsBeginIter;
      ++degreesIter;
    }
  }

//  int getKnotIdx(int dim, int degree, T u) const {
//    auto idx = degree + (int)(u / knotDeltas_[dim]);
//    if (idx == pointDims_[dim]) {
//      idx = pointDims_[dim] - 1;
//    }
//    assert(idx >= degree && idx < pointDims_[dim]);
//    return idx;
//  }

  int getKnotIdx(int dim, int degree, double u) const {
    auto idx = (degree) + (int)(u / knotDeltas_[dim]);
    if (idx == (pointDims_[dim])) {
      idx = (pointDims_[dim] - 1);
    }
    assert(idx >= (degree) && idx < (pointDims_[dim]));
    return idx;
  }

//  template <typename OutIter>
//  void computeBasis(int k, int dim, int degree, T u, OutIter n) const {
//    n[0] = T(1.0);
//
//    if (k < 2 * degree - 1 || k >= pointDims_[dim] - degree + 1) {
//      //this can be done more efficient
//      const auto& knot = knots_[dim];
//      for (int i = 1; i <= degree; ++i) {
//        n[i] = T(0.0);
//
//        for (int j = i - 1; j >= 0; --j) {
//          T a = T(0.0);
//          if (knot[k + i - j] != knot[k - j]) {
//            a = (u - knot[k - j]) / (knot[k + i - j] - knot[k - j]);
//          }
//
//          n[j + 1] += n[j] * (T(1.0) - a);
//          n[j] = n[j] * a;
//        }
//      }
//    } else {
//      auto m = u / knotDeltas_[dim] - (T(k - degree));
//
//      for (int i = 1; i <= degree; ++i) {
//        n[i] = T(0.0);
//
//        for (int j = i - 1; j >= 0; --j) {
//          auto a = (m + T(j)) / T(i);
//          n[j + 1] += n[j] * (T(1.0) - a);
//          n[j] = n[j] * a;
//        }
//      }
//    }
//  }

  template <typename OutIter>
    void computeBasis(int k, int dim, int degree, double u, OutIter n) const {
      n[0] = T(1.0);

      if (k < 2 * degree - 1 || k >= pointDims_[dim] - degree + 1) {
        //this can be done more efficient
        const auto& knot = knots_[dim];
        for (int i = 1; i <= degree; ++i) {
          n[i] = T(0.0);

          for (int j = i - 1; j >= 0; --j) {
            T a = T(0.0);
            if (knot[k + i - j] != knot[k - j]) {
              a = (u - knot[k - j]) / (knot[k + i - j] - knot[k - j]);
            }

            n[j + 1] += n[j] * (T(1.0) - a);
            n[j] = n[j] * a;
          }
        }
      } else {
        auto m = u / knotDeltas_[dim] - (T(k - degree));

        for (int i = 1; i <= degree; ++i) {
          n[i] = T(0.0);

          for (int j = i - 1; j >= 0; --j) {
            auto a = (m + T(j)) / T(i);
            n[j + 1] += n[j] * (T(1.0) - a);
            n[j] = n[j] * a;
          }
        }
      }
    }
//  template <typename OutIter>
//  void computeBasis(T k, int dim, int degree, T u, OutIter n) const {
//    n[0] = T(1.0);
//
//    if (k < T(2 * degree - 1) || k >= T(pointDims_[dim] - degree + 1)) {
//      //this can be done more efficient
//      const auto& knot = knots_[dim];
//      for (T i = T(1); i <= T(degree); i + T(1)) {
//        n[i] = T(0.0);
//
//        for (T j = i - T(1); j >= T(0); j - T(1)) {
//          T a = T(0.0);
//          if (knot[k + i - j] != knot[k - j]) {
//            a = (u - knot[k - j]) / (knot[k + i - j] - knot[k - j]);
//          }
//
//          n[j + 1] += n[j] * (T(1.0) - a);
//          n[j] = n[j] * a;
//        }
//      }
//    } else {
//      auto m = u / knotDeltas_[dim] - (T(k - degree));
//
//      for (T i = T(1); i <= T(degree); i + T(1)) {
//        n[i] = T(0.0);
//
//        for (T j = i - T(1); j >= 0; j - T(1)) {
//          auto a = (m + T(j)) / T(i);
//          n[j + 1] += n[j] * (T(1.0) - a);
//          n[j] = n[j] * a;
//        }
//      }
//    }
//  }


  T getMinU(int /*inputD*/) const {
    return T(0.0);
  }

  T getMaxU(int /*inputD*/) const {
    return T(1.0);
  }

  template <typename Archive>
  void serialize(Archive& archive) {
    archive(pointDims_, degrees_, knotDeltas_, knots_);
  }

private:
  std::array<int, inputDim> pointDims_;
  std::array<int, inputDim> degrees_;
  std::array<double, inputDim> knotDeltas_;

  KnotsT knots_;
};

template <typename T, std::size_t inputDim, std::size_t outputDim, WeightType weightType>
class NurbsWeights {
};

template <typename T, std::size_t inputDim, std::size_t outputDim>
class NurbsWeights < T, inputDim, outputDim, WeightType::RATIONAL > {
private:
  typedef boost::multi_array<T, inputDim> WeightsT;

public:
  NurbsWeights() = default;

  NurbsWeights(const NurbsWeights<T, inputDim, outputDim, WeightType::RATIONAL>& o) {
    assignMultiArray(o.getWeights(), weights_);
  }
  NurbsWeights(const NurbsWeights<T, inputDim, outputDim, WeightType::NON_RATIONAL>& o) {
    assignMultiArray(o.getWeights(), weights_);
  }

  NurbsWeights<T, inputDim, outputDim, WeightType::RATIONAL> operator=(const NurbsWeights<T, inputDim, outputDim, WeightType::RATIONAL>& o) {
    if (&o != this) {
      assignMultiArray(o.getWeights(), weights_);
    }
    return *this;
  }
  NurbsWeights<T, inputDim, outputDim, WeightType::RATIONAL> operator=(const NurbsWeights<T, inputDim, outputDim, WeightType::NON_RATIONAL>& o) {
    assignMultiArray(o.getWeights(), weights_);
    return *this;
  }

  const WeightsT& getWeights() const {
    return weights_;
  }

  WeightsT& getWeights() {
    return weights_;
  }

  void create(const WeightsT& weights) {
    assignMultiArray(weights, weights_);
  }

  template <typename Iter>
  void create(Iter numPointsBeginIter) {
    std::array<int, inputDim> numPoints;
    std::copy_n(numPointsBeginIter, inputDim, numPoints.begin());

    weights_.resize(numPoints);
    std::fill(weights_.data(), weights_.data() + weights_.num_elements(), T(1.0));
  }

  template <typename NsT, typename KsT, typename DegreesT, typename PointsT, typename OutIter>
  void getPoint(const NsT& ns, const KsT& ks, const DegreesT& degrees, const PointsT& points, OutIter pBegin) const {
    T weightSum = T(0.0);
    std::fill(pBegin, pBegin + outputDim, T(0.0));

    applyGetPoint(ns, ks, degrees, [&, this] (const std::array<int, inputDim>& idx, T val) {
      val *= weights_(idx);
      weightSum += val;

      auto resP = pBegin;
      for (int i = 0; i < (int)outputDim; ++i) {
        *resP += val * points(idx)[i];
        ++resP;
      }
    });

    auto resP = pBegin;

                if(weightSum == T(0)) weightSum = T(1);

    for (int i = 0; i < (int)outputDim; ++i) {
      *resP /= weightSum;
      ++resP;
    }
  }

  template <typename Archive>
  void serialize(Archive& archive) {
    archive(weights_);
  }

private:
  WeightsT weights_;
};

template <typename T, std::size_t inputDim, std::size_t outputDim>
class NurbsWeights < T, inputDim, outputDim, WeightType::NON_RATIONAL > {
private:
  typedef boost::multi_array<T, inputDim> WeightsT;

public:
  WeightsT getWeights() const {
    WeightsT weights;
    weights.resize(numPoints_);
    std::fill(weights.data(), weights.data() + weights.num_elements(), T(1.0));
    return weights;
  }

  template <typename Iter>
  void create(Iter numPointsBeginIter) {
    auto iter = numPoints_.begin();
    for (std::size_t i = 0; i < inputDim; ++i, ++iter, ++numPointsBeginIter) {
      *iter = static_cast<int>(*numPointsBeginIter);
    }
  }

  template <typename NsT, typename KsT, typename DegreesT, typename PointsT, typename OutIter>
  void getPoint(const NsT& ns, const KsT& ks, const DegreesT& degrees, const PointsT& points, OutIter pBegin) const {
    std::fill(pBegin, pBegin + outputDim, T(0.0));

    applyGetPoint(ns, ks, degrees, [&, this] (const std::array<int, inputDim>& idx, T val) {
      auto resP = pBegin;
      for (int i = 0; i < (int)outputDim; ++i) {
        *resP += val * points(idx)[i];
        ++resP;
      }
    });
  }

  template <typename Archive>
  void serialize(Archive& archive) {
    archive(numPoints_);
  }

private:
  std::array<int, inputDim> numPoints_;
};

template <typename T, std::size_t inputDim, std::size_t outputDim, KnotType knotType, WeightType weightType>
class Nurbs {
  typedef std::array<int, inputDim> DegreesT;
  typedef boost::multi_array<std::array<T, outputDim>, inputDim> PointsT;

  typedef NurbsKnots<T, inputDim, outputDim, knotType> NurbsKnotsT;
  typedef NurbsWeights<T, inputDim, outputDim, weightType> NurbsWeightsT;

  typedef std::array<std::vector<T>, inputDim> NsT;
  typedef std::array<int, inputDim> KsT;

  NurbsKnotsT knots_;
  NurbsWeightsT weights_;

  DegreesT degrees_; //one degree for each input dimension
  PointsT points_; //one n-dimensional point for each grid point

public:
  template <std::size_t inputDim1, WeightType weightType1>
  friend struct NurbsCeresPrivate::CreateParameterVector;

  template <std::size_t inputDim1, std::size_t outputDim1, KnotType knotType1, WeightType weightType1>
  friend class NurbsCeres;

  template <typename T1, std::size_t inputDim1, std::size_t outputDim1, KnotType knotType1, WeightType weightType1>
  friend class Nurbs;

  Nurbs() = default;

  Nurbs(const Nurbs<T, inputDim, outputDim, knotType, weightType>& o) : knots_(o.knots_), weights_(o.weights_), degrees_(o.degrees_) {
    assignMultiArray(o.points_, points_);
      assert(validDimensions() == true);
  }

  template <KnotType knotType2, WeightType weightType2>
  Nurbs(const Nurbs<T, inputDim, outputDim, knotType2, weightType2>& o) : knots_(o.knots_), weights_(o.weights_), degrees_(o.degrees_) {
    assignMultiArray(o.points_, points_);
      assert(validDimensions() == true);
  }

  Nurbs<T, inputDim, outputDim, knotType, weightType>& operator=(const Nurbs<T, inputDim, outputDim, knotType, weightType>& o) {
    if (&o != this) {
      knots_ = o.knots_;
      weights_ = o.weights_;
      degrees_ = o.degrees_;
      assignMultiArray(o.points_, points_);
        assert(validDimensions() == true);
    }
    return *this;
  }
  NurbsWeightsT& getWeights()
  {
    return weights_;
  }

  NurbsWeightsT& getWeights() const
  {
    return weights_;
  }

  template <typename T2, std::size_t inputDim2, std::size_t outputDim2, KnotType knotType2, WeightType weightType2>
  Nurbs<T, inputDim, outputDim, knotType, weightType>& operator=(const Nurbs<T2, inputDim2, outputDim2, knotType2, weightType2>& o) {
    knots_ = o.knots_;
    weights_ = o.weights_;
    degrees_ = o.degrees_;
    assignMultiArray(o.points_, points_);
      assert(validDimensions() == true);

    return *this;
  }

  const PointsT& getPoints() const {
    return points_;
  }

  PointsT& getPoints() {
    return points_;
  }

  void init(int degree, const NurbsKnotsT& knots, const NurbsWeightsT& weights, const PointsT& points) {
    DegreesT degrees;
    std::fill(degrees.begin(), degrees.end(), degree);
    init(degrees, knots, weights, points);
  }

  void init(const DegreesT& degrees, const NurbsKnotsT& knots, const NurbsWeightsT& weights, const PointsT& points) {
    degrees_ = degrees;
    assignMultiArray(points, points_);
    knots_ = knots;
    weights_ = weights;

     assert(validDimensions() == true);
  }

  T getMinU(int inputD) const {
    return knots_.getMinU(inputD);
  }

  T getMaxU(int inputD) const {
    return knots_.getMaxU(inputD);
  }

  template <typename InIter, typename OutIter>
  void getPoint(InIter uBegin, OutIter rBegin) const {
    KsT ks;
    NsT ns;

    getBasis(uBegin, ks.begin(), ns.begin());
    weights_.getPoint(ns, ks, degrees_, points_, rBegin);
  }

  template <typename InIter, typename IdxIter, typename BaseIter>
  void getBasis(InIter uBegin, IdxIter idxIter, BaseIter baseIter) const {
    for (int dim = 0; dim < (int)inputDim; ++dim) {
      const auto& degree = degrees_[dim];
      const auto& u = *uBegin;

      auto& k = *idxIter;
      auto& n = *baseIter;

      n.resize(degree + 1);

      k = knots_.getKnotIdx(dim, degree, u);
      //T k_templated = knots_.getKnotIdx(dim, degree, u);;
      knots_.computeBasis(k, dim, degree, u, n.begin());
      assert(isBaseCorrect(n.begin(), n.end()));

      ++uBegin;
      ++idxIter;
      ++baseIter;
    }
  }

  bool validDimensions() const {
    auto pointsShape = points_.shape();
    auto weights = weights_.getWeights();
    auto weightsShape = weights.shape();
    auto knots = knots_.getKnot();

    for (int dim = 0; dim < (int)inputDim; ++dim) {
      if ((int)knots[dim].size() < 2 * degrees_[dim]) {
        return false;
      }

      if ((int)knots[dim].size() != (int)pointsShape[dim] + degrees_[dim] + 1) {
        return false;
      }

      if ((int)knots[dim].size() != (int)weightsShape[dim] + degrees_[dim] + 1) {
        return false;
      }
    }

    return true;
  }

  bool validKnotVector() const {
    auto knots = knots_.getKnots();

    for (int dim = 0; dim < (int)inputDim; ++dim) {
      if (std::is_sorted(knots[dim].begin(), knots[dim].end()) == false) {
        return false;
      }
    }

    return true;
  }

  template <typename Archive>
  void serialize(Archive& archive) {
    archive(knots_, weights_, degrees_, points_);
  }
};
}
}




#endif /* CARTOGRAPHER_MAPPING_3D_SCAN_MATCHING_NURBS_H_ */
