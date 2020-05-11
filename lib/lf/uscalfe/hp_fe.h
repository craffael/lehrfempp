#ifndef LF_USCALFE_HP_FE_H_
#define LF_USCALFE_HP_FE_H_

/**
 * @file
 * @brief Data structures representing HP finite elements
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include "lagr_fe.h"

#include <cmath>
#include <memory>
#include <vector>


namespace lf::uscalfe {
/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::assemble::ldof_idx_t;
/** Type for vector length/matrix sizes */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;


/**
 * @brief Struct for computing Legrendre and integrated Legendre polynomials
 */
template<typename SCALAR>
struct LegendrePoly {
    static SCALAR eval(unsigned n, double x) {
	double Ljm1 = 1;
	double Lj = x;
	if (n == 0) {
	    return Ljm1;
	}
	else if (n == 1) {
	    return Lj;
	}
	else {
	    for (unsigned j = 1 ; j < n ; ++j) {
		double Ljp1 = ((2*j+1)*x*Lj - j*Ljm1) / (j+1);
		Ljm1 = Lj;
		Lj = Ljp1;
	    }
	    return Lj;
	}
    }

    static SCALAR integral(unsigned n, double x) {
	if (n == 0) {
	    return -1;
	}
	else if (n == 1) {
	    return x;
	}
	else {
	    double Ljm2 = 1;
	    double Ljm1 = x;
	    double Lj = (3*x*x - 1) / 2;
	    for (unsigned j = 2 ; j < n ; ++j) {
		double Ljp1 = ((2*j+1)*x*Lj - j*Ljm1) / (j+1);
		Ljm2 = Ljm1;
		Ljm1 = Lj;
		Lj = Ljp1;
	    }
	    return (Lj - Ljm2) / (2*n-1);
	}
    }
};


/**
 * @brief Computes Chebyshev interpolation nodes in [0, 1]
 * @param n Degree of the Chebyshev interpolation nodes
 * @returns An Eigen vector containing the interpolation nodes on [0, 1]
 */
Eigen::VectorXd chebyshevNodes(unsigned n);



/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief HP Finite Elements of arbitrary degree on segments
 *
 * The Shape Functions are taken from the following paper: https://arxiv.org/pdf/1504.03025.pdf
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeHPSegment final : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHPSegment(const FeHPSegment &) = default;
  FeHPSegment(FeHPSegment &&) = default;
  FeHPSegment &operator=(const FeHPSegment &) = default;
  FeHPSegment &operator=(FeHPSegment &&) = default;
  ~FeHPSegment() = default;

  FeHPSegment(unsigned degree, nonstd::span<const lf::mesh::Orientation> rel_orient) : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(), degree_(degree), rel_orient_(rel_orient) {
    eval_nodes_ = ComputeEvaluationNodes();
  }

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return degree_; }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return degree_ + 1;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    return codim == 0 ? degree_ - 1 : 1;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,
   * sub_idx_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim,
      sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
      LF_ASSERT_MSG(refcoords.rows() == 1, "refcoords must be a row vector");
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), refcoords.cols());
      // Get the shape functions associated with the vertices
      result.row(0) = refcoords.unaryExpr([&](double x) -> SCALAR { return 1-x; });
      result.row(1) = refcoords.unaryExpr([&](double x) -> SCALAR { return x; });
      // Get the shape functions associated with the interior of the segment
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  result.row(i+2) = refcoords.unaryExpr([&](double x) -> SCALAR {
	    return LegendrePoly<SCALAR>::integral(i+2, 1-2*x);
	  });
      }
      return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
      LF_ASSERT_MSG(refcoords.rows() == 1, "refcoords must be a row vector");
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), refcoords.cols());
      // Get the gradient of the shape functions associated with the vertices
      result.row(0) = refcoords.unaryExpr([&](double x) -> SCALAR { return -1; });
      result.row(1) = refcoords.unaryExpr([&](double x) -> SCALAR { return 1; });
      // Get the shape functions associated with the interior of the segment
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  result.row(i+2) = refcoords.unaryExpr([&](double x) -> SCALAR {
	    return -2*LegendrePoly<SCALAR>::eval(i+1, 1-2*x);
	  });
      }
      return result;
  }

  /**
   * @brief Evaluation nodes are the endpoints of the segment and the Chebyshev
   * nodes of degree p-1 on the segment
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
      return eval_nodes_;
  }

  /**
   * @brief p+1 shape functions
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
      return degree_ + 1;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodevals) const override {
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> shape_functions_at_nodes = EvalReferenceShapeFunctions(EvaluationNodes());
    return shape_functions_at_nodes.transpose().fullPivHouseholderQr().solve(nodevals.transpose()).transpose();
  }

 private:
  unsigned degree_;
  Eigen::MatrixXd eval_nodes_;
  nonstd::span<const lf::mesh::Orientation> rel_orient_;

  Eigen::MatrixXd ComputeEvaluationNodes() const {
      Eigen::MatrixXd nodes(1, degree_+1);
      nodes(0, 0) = 0;
      nodes(0, 1) = 1;
      if (degree_ > 1) {
	  nodes.block(0, 2, 1, degree_-1) = chebyshevNodes(degree_-1).transpose();
      }
      return nodes;
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief HP Finite Elements of arbitrary degree on triangles
 *
 * The Shape Functions are taken from the following paper: https://arxiv.org/pdf/1504.03025.pdf
 * where Legendre polynomials are taken for the face bubbles instead of Jacobi polynomials
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeHPTria final : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHPTria(const FeHPTria &) = default;
  FeHPTria(FeHPTria &&) = default;
  FeHPTria &operator=(const FeHPTria &) = default;
  FeHPTria &operator=(FeHPTria &&) = default;
  ~FeHPTria() = default;

  FeHPTria(unsigned degree, nonstd::span<const lf::mesh::Orientation> rel_orient)
      : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(), degree_(degree), rel_orient_(rel_orient), eval_nodes_() {
    eval_nodes_ = ComputeEvaluationNodes();
  }

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return degree_; }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return (degree_ + 1) * (degree_ + 2) / 2;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    switch (codim) {
      case 0:
        if (degree_ <= 2) {
          return 0;
        } else {
          return (degree_ - 2) * (degree_ - 1) / 2;
        }
      case 1:
        return degree_ - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
    }
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim,
      sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), refcoords.cols());
      // Compute the barycentric coordinate functions
      const Eigen::RowVectorXd l1 = Eigen::RowVectorXd::Ones(refcoords.cols()) - refcoords.row(0) - refcoords.row(1);
      const Eigen::RowVectorXd l2 = refcoords.row(0);
      const Eigen::RowVectorXd l3 = refcoords.row(1);
      // Get the basis functions associated with the vertices
      result.row(0) = l1.unaryExpr([&](double x) -> SCALAR { return x; });
      result.row(1) = l2.unaryExpr([&](double x) -> SCALAR { return x; });
      result.row(2) = l3.unaryExpr([&](double x) -> SCALAR { return x; });
      // Get the basis functions associated with the first edge
      Eigen::RowVectorXd l1n(refcoords.cols());
      for (int i = 0 ; i < refcoords.cols() ; ++i) {
	  l1n[i] = l1[i]+l2[i] == 0 ? SCALAR(0) : l1[i]/(l1[i]+l2[i]);
      }
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  if (rel_orient_[0] == lf::mesh::Orientation::positive) {
	      result.row(3+i) = ((l1 + l2).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
				 l1n.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); })).matrix();
	  }
	  else {
	      result.row(degree_+1-i) = ((l1 + l2).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
				         l1n.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 1-2*x); })).matrix();
	  }
      }
      // Get the basis functions associated with the second edge
      Eigen::RowVectorXd l2n(refcoords.cols());
      for (int i = 0 ; i < refcoords.cols() ; ++i) {
	  l2n[i] = l2[i]+l3[i] == 0 ? SCALAR(0) : l2[i]/(l2[i]+l3[i]);
      }
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  if (rel_orient_[1] == lf::mesh::Orientation::positive) {
	      result.row(degree_+2+i) = ((l2 + l3).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
					 l2n.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); })).matrix();
	  }
	  else {
	      result.row(2*degree_-i) = ((l2 + l3).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
					 l2n.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 1-2*x); })).matrix();
	  }
      }
      // Get the basis functions associated with the third edge
      Eigen::RowVectorXd l3n(refcoords.cols());
      for (int i = 0 ; i < refcoords.cols() ; ++i) {
	  l3n[i] = l3[i]+l1[i] == 0 ? SCALAR(0) : l3[i]/(l3[i]+l1[i]);
      }
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  if (rel_orient_[2] == lf::mesh::Orientation::positive) {
	      result.row(2*degree_+1+i) = ((l3 + l1).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
					   l3n.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); })).matrix();
	  }
	  else {
	      result.row(3*degree_-1-i) = ((l3 + l1).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
					   l3n.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 1-2*x); })).matrix();
	  }
      }
      // Get the basis functions associated with the interior of the triangle
      if (degree_ > 2) {
	  int idx = 3 * degree_;
	  for (int i = 0 ; i < degree_-2 ; ++i) {
	      for (int j = 0 ; j < degree_-i-2 ; ++j) {
		  result.row(idx) = ((l2 + l3).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
				     (l2.array()/(l2+l3).array()).unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); }) *
				     l1.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(j+2, 2*x-1); })).matrix();
		  ++ idx;
	      }
	  }
      }
      return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), 2 * refcoords.cols());
    // Compute the barycentric coordinate functions
    const Eigen::RowVectorXd l1 = Eigen::RowVectorXd::Ones(refcoords.cols()) - refcoords.row(0) - refcoords.row(1);
    const Eigen::RowVectorXd l2 = refcoords.row(0);
    const Eigen::RowVectorXd l3 = refcoords.row(1);
    const Eigen::RowVectorXd l1_dx = Eigen::RowVectorXd::Constant(refcoords.cols(), -1);
    const Eigen::RowVectorXd l1_dy = Eigen::RowVectorXd::Constant(refcoords.cols(), -1);
    const Eigen::RowVectorXd l2_dx = Eigen::RowVectorXd::Constant(refcoords.cols(), 1);
    const Eigen::RowVectorXd l2_dy = Eigen::RowVectorXd::Constant(refcoords.cols(), 0);
    const Eigen::RowVectorXd l3_dx = Eigen::RowVectorXd::Constant(refcoords.cols(), 0);
    const Eigen::RowVectorXd l3_dy = Eigen::RowVectorXd::Constant(refcoords.cols(), 1);
    for (int i = 0 ; i < refcoords.cols() ; ++i) {
      // Get the gradient of the basis functions associated with the vertices
      result(0, 2*i+0) = l1_dx[i];
      result(0, 2*i+1) = l1_dy[i];
      result(1, 2*i+0) = l2_dx[i];
      result(1, 2*i+1) = l2_dy[i];
      result(2, 2*i+0) = l3_dx[i];
      result(2, 2*i+1) = l3_dy[i];
      // Get the gradient of the basis functions associated with the first edge
      const SCALAR l1p2 = l1[i] + l2[i];
      const SCALAR l1p2_dx = l1_dx[i] + l2_dx[i];
      const SCALAR l1p2_dy = l1_dy[i] + l2_dy[i];
      const SCALAR l1n = l1p2 == 0 ? SCALAR(0) : (l1[i] / l1p2);
      const SCALAR l1n_dx = l1p2 == 0 ? SCALAR(0) : ((l1_dx[i]*l1p2 - l1[i]*l1p2_dx) / (l1p2 * l1p2));
      const SCALAR l1n_dy = l1p2 == 0 ? SCALAR(0) : ((l1_dy[i]*l1p2 - l1[i]*l1p2_dy) / (l1p2 * l1p2));
      for (int j = 0 ; j < degree_-1 ; ++j) {
	  if (rel_orient_[0] == lf::mesh::Orientation::positive) {
	      const SCALAR leg1inte = LegendrePoly<SCALAR>::integral(j+2, 2*l1n-1);
	      const SCALAR leg1eval = LegendrePoly<SCALAR>::eval(j+1, 2*l1n-1);
	      result(3+j, 2*i+0) = l1p2_dx*(j+2)*std::pow(l1p2, j+1)*leg1inte + std::pow(l1p2, j+2)*2*l1n_dx*leg1eval;
	      result(3+j, 2*i+1) = l1p2_dy*(j+2)*std::pow(l1p2, j+1)*leg1inte + std::pow(l1p2, j+2)*2*l1n_dy*leg1eval;
	  }
	  else {
	      const SCALAR leg1inte = LegendrePoly<SCALAR>::integral(j+2, 1-2*l1n);
	      const SCALAR leg1eval = LegendrePoly<SCALAR>::eval(j+1, 1-2*l1n);
	      result(degree_+1-j, 2*i+0) = l1p2_dx*(j+2)*std::pow(l1p2, j+1)*leg1inte - std::pow(l1p2, j+2)*2*l1n_dx*leg1eval;
	      result(degree_+1-j, 2*i+1) = l1p2_dy*(j+2)*std::pow(l1p2, j+1)*leg1inte - std::pow(l1p2, j+2)*2*l1n_dy*leg1eval;
	  }
      }
      // Get the gradient of the basis functions associated with the second edge
      const SCALAR l2p3 = l2[i] + l3[i];
      const SCALAR l2p3_dx = l2_dx[i] + l3_dx[i];
      const SCALAR l2p3_dy = l2_dy[i] + l3_dy[i];
      const SCALAR l2n = l2p3 == 0 ? SCALAR(0) : (l2[i] / l2p3);
      const SCALAR l2n_dx = l2p3 == 0 ? SCALAR(0) : ((l2_dx[i]*l2p3 - l2[i]*l2p3_dx) / (l2p3 * l2p3));
      const SCALAR l2n_dy = l2p3 == 0 ? SCALAR(0) : ((l2_dy[i]*l2p3 - l2[i]*l2p3_dy) / (l2p3 * l2p3));
      for (int j = 0 ; j < degree_-1 ; ++j) {
	  if (rel_orient_[1] == lf::mesh::Orientation::positive) {
	      const SCALAR leg2inte = LegendrePoly<SCALAR>::integral(j+2, 2*l2n-1);
	      const SCALAR leg2eval = LegendrePoly<SCALAR>::eval(j+1, 2*l2n-1);
	      result(2+degree_+j, 2*i+0) = l2p3_dx*(j+2)*std::pow(l2p3, j+1)*leg2inte + std::pow(l2p3, j+2)*2*l2n_dx*leg2eval;
	      result(2+degree_+j, 2*i+1) = l2p3_dy*(j+2)*std::pow(l2p3, j+1)*leg2inte + std::pow(l2p3, j+2)*2*l2n_dy*leg2eval;
	  }
	  else {
	      const SCALAR leg2inte = LegendrePoly<SCALAR>::integral(j+2, 1-2*l2n);
	      const SCALAR leg2eval = LegendrePoly<SCALAR>::eval(j+1, 1-2*l2n);
	      result(2*degree_-j, 2*i+0) = l2p3_dx*(j+2)*std::pow(l2p3, j+1)*leg2inte - std::pow(l2p3, j+2)*2*l2n_dx*leg2eval;
	      result(2*degree_-j, 2*i+1) = l2p3_dy*(j+2)*std::pow(l2p3, j+1)*leg2inte - std::pow(l2p3, j+2)*2*l2n_dy*leg2eval;
	  }
      }
      // Get the gradient of the basis functions associated with the third edge
      const SCALAR l3p1 = l3[i] + l1[i];
      const SCALAR l3p1_dx = l3_dx[i] + l1_dx[i];
      const SCALAR l3p1_dy = l3_dy[i] + l1_dy[i];
      const SCALAR l3n = l3p1 == 0 ? SCALAR(0) : (l3[i] / l3p1);
      const SCALAR l3n_dx = l3p1 == 0 ? SCALAR(0) : ((l3_dx[i]*l3p1 - l3[i]*l3p1_dx) / (l3p1 * l3p1));
      const SCALAR l3n_dy = l3p1 == 0 ? SCALAR(0) : ((l3_dy[i]*l3p1 - l3[i]*l3p1_dy) / (l3p1 * l3p1));
      for (int j = 0 ; j < degree_-1 ; ++j) {
	  if (rel_orient_[2] == lf::mesh::Orientation::positive) {
	      const SCALAR leg3inte = LegendrePoly<SCALAR>::integral(j+2, 2*l3n-1);
	      const SCALAR leg3eval = LegendrePoly<SCALAR>::eval(j+1, 2*l3n-1);
	      result(1+2*degree_+j, 2*i+0) = l3p1_dx*(j+2)*std::pow(l3p1, j+1)*leg3inte + std::pow(l3p1, j+2)*2*l3n_dx*leg3eval;
	      result(1+2*degree_+j, 2*i+1) = l3p1_dy*(j+2)*std::pow(l3p1, j+1)*leg3inte + std::pow(l3p1, j+2)*2*l3n_dy*leg3eval;
	  }
	  else {
	      const SCALAR leg3inte = LegendrePoly<SCALAR>::integral(j+2, 1-2*l3n);
	      const SCALAR leg3eval = LegendrePoly<SCALAR>::eval(j+1, 1-2*l3n);
	      result(3*degree_-1-j, 2*i+0) = l3p1_dx*(j+2)*std::pow(l3p1, j+1)*leg3inte - std::pow(l3p1, j+2)*2*l3n_dx*leg3eval;
	      result(3*degree_-1-j, 2*i+1) = l3p1_dy*(j+2)*std::pow(l3p1, j+1)*leg3inte - std::pow(l3p1, j+2)*2*l3n_dy*leg3eval;
	  }
      }
      // Get the gradient of the basis functions associated with the interior of the triangle
      if (degree_ > 2) {
	  int idx = 3 * degree_;
	  for (int j = 0 ; j < degree_-2 ; ++j) {
	      SCALAR legjinte = LegendrePoly<SCALAR>::integral(j+2, 2*l2n-1);
	      SCALAR legjeval = LegendrePoly<SCALAR>::eval(j+1, 2*l2n-1);
	      for (int k = 0 ; k < degree_-j-2 ; ++k) {
		  SCALAR legkinte = LegendrePoly<SCALAR>::integral(k+2, 2*l1[i]-1);
		  SCALAR legkeval = LegendrePoly<SCALAR>::eval(k+1, 2*l1[i]-1);
		  result(idx, 2*i+0) = (j+2)*std::pow(l2p3, j+1)*legjinte*legkinte + 2*l3[i]*std::pow(l2p3, j)*legjeval*legkinte - 2*std::pow(l2p3, j+2)*legjinte*legkeval;
		  result(idx, 2*i+1) = (j+2)*std::pow(l2p3, j+1)*legjinte*legkinte - 2*l2[i]*std::pow(l2p3, j)*legjeval*legkinte - 2*std::pow(l2p3, j+2)*legjinte*legkeval;
		  ++idx;
	      }
	  }
      }
    }
    return result;
  }

  /**
   * @brief Evaluation nodes are the vertices, the Chebyshev nodes of degree p-1
   * on the edges and the corresponding nodes on the triangle
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return eval_nodes_;
  }

  /**
   * @brief (p+1)*(p+2)/2 evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodevals) const override {
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> shape_functions_at_nodes = EvalReferenceShapeFunctions(EvaluationNodes());
    return shape_functions_at_nodes.transpose().fullPivHouseholderQr().solve(nodevals.transpose()).transpose();
  }

 private:
  unsigned degree_;
  Eigen::MatrixXd eval_nodes_;
  nonstd::span<const lf::mesh::Orientation> rel_orient_;

  Eigen::MatrixXd ComputeEvaluationNodes() const {
    Eigen::MatrixXd eval_nodes(2, (degree_ + 1) * (degree_ + 2) / 2);
    const auto cheb = chebyshevNodes(degree_-1);
    // Add the evaluation nodes corresponding to the vertices of the triangle
    eval_nodes(0, 0) = 0;
    eval_nodes(1, 0) = 0;
    eval_nodes(0, 1) = 1;
    eval_nodes(1, 1) = 0;
    eval_nodes(0, 2) = 0;
    eval_nodes(1, 2) = 1;
    // Add the evaluation nodes corresponding to the edges of the triangle
    for (int i = 0; i < degree_ - 1; ++i) {
      eval_nodes(0, 3 + i) = cheb[i];
      eval_nodes(1, 3 + i) = 0;
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      eval_nodes(0, 2 + degree_ + i) = 1. - cheb[i];
      eval_nodes(1, 2 + degree_ + i) = cheb[i];
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      eval_nodes(0, 1 + 2*degree_ + i) = 0;
      eval_nodes(1, 1 + 2*degree_ + i) = 1. - cheb[i];
    }
    // Add the evaluation nodes corresponding to the interior of the triangle
    if (degree_ > 2) {
      int idx = 3 * degree_;
      for (int i = 0; i < degree_ - 2; ++i) {
        for (int j = 0; j < degree_-2-i ; ++j) {
          eval_nodes(0, idx) = cheb[j];
          eval_nodes(1, idx) = cheb[i];
	  ++idx;
        }
      }
    }
    return eval_nodes;
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief HP Finite Elements of arbitrary degre
 * e on quadrilaterals
 *
 * The Shape Functions are taken from the following paper: https://arxiv.org/pdf/1504.03025.pdf
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeHPQuad final : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHPQuad(const FeHPQuad &) = default;
  FeHPQuad(FeHPQuad &&) = default;
  FeHPQuad &operator=(const FeHPQuad &) = default;
  FeHPQuad &operator=(FeHPQuad &&) = default;
  ~FeHPQuad() = default;

  FeHPQuad(unsigned degree, nonstd::span<const lf::mesh::Orientation> rel_orient)
      : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(), degree_(degree), rel_orient_(rel_orient), eval_nodes_(), fe1d_(degree, rel_orient) {
    eval_nodes_ = ComputeEvaluationNodes();
  }

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return degree_; }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return (degree_ + 1) * (degree_ + 1);
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    switch (codim) {
      case 0:
        return (degree_ - 1) * (degree_ - 1);
      case 1:
        return degree_ - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
    }
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim,
      sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), refcoords.cols());
    // Compute the 1D shape functions at the x and y coordinates
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_x = fe1d_.EvalReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_y = fe1d_.EvalReferenceShapeFunctions(refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_x = fe1d_.EvalReferenceShapeFunctions(Eigen::RowVectorXd::Constant(refcoords.cols(), 1)-refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_y = fe1d_.EvalReferenceShapeFunctions(Eigen::RowVectorXd::Constant(refcoords.cols(), 1)-refcoords.row(1));
    // Get the basis functions associated with the vertices
    result.row(0) = (sf1d_x.row(0).array() * sf1d_y.row(0).array()).matrix();
    result.row(1) = (sf1d_x.row(1).array() * sf1d_y.row(0).array()).matrix();
    result.row(2) = (sf1d_x.row(1).array() * sf1d_y.row(1).array()).matrix();
    result.row(3) = (sf1d_x.row(0).array() * sf1d_y.row(1).array()).matrix();
    // Get the basis functions associated with the first edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	if (rel_orient_[0] == lf::mesh::Orientation::positive) {
	    result.row(4+i) = (sf1d_x.row(2+i).array() * sf1d_y.row(0).array()).matrix();
	}
	else {
	    result.row(2+degree_-i) = (sf1df_x.row(2+i).array() * sf1d_y.row(0).array()).matrix();
	}
    }
    // Get the basis functions associated with the second edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	if (rel_orient_[1] == lf::mesh::Orientation::positive) {
	    result.row(3+degree_+i) = (sf1d_x.row(1).array() * sf1d_y.row(2+i).array()).matrix();
	}
	else {
	    result.row(1+2*degree_-i) = (sf1d_x.row(1).array() * sf1df_y.row(2+i).array()).matrix();
	}
    }
    // Get the basis functions associated with the third edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	if (rel_orient_[2] == lf::mesh::Orientation::positive) {
	    result.row(2+2*degree_+i) = (sf1df_x.row(2+i).array() * sf1d_y.row(1).array()).matrix();
	}
	else {
	    result.row(3*degree_-i) = (sf1d_x.row(2+i).array() * sf1d_y.row(1).array()).matrix();
	}
    }
    // Get the basis functions associated with the fourth edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	if (rel_orient_[3] == lf::mesh::Orientation::positive) {
	    result.row(1+3*degree_+i) = (sf1d_x.row(0).array() * sf1df_y.row(2+i).array()).matrix();
	}
	else {
	    result.row(4*degree_-1-i) = (sf1d_x.row(0).array() * sf1d_y.row(2+i).array()).matrix();
	}
    }
    // Get the basis functions associated with the interior of the quad
    for (int i = 0 ; i < degree_-1 ; ++i) {
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    result.row(4*degree_ + (degree_-1)*i + j) = (sf1d_x.row(j+2).array() * sf1d_y.row(i+2).array()).matrix();
	}
    }
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), 2*refcoords.cols());
    // Compute the gradient of the 1D shape functions at the x and y coordinates
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_x = fe1d_.EvalReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_y = fe1d_.EvalReferenceShapeFunctions(refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_dx = fe1d_.GradientsReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_dy = fe1d_.GradientsReferenceShapeFunctions(refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_x = fe1d_.EvalReferenceShapeFunctions(Eigen::RowVectorXd::Constant(refcoords.cols(), 1) - refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_y = fe1d_.EvalReferenceShapeFunctions(Eigen::RowVectorXd::Constant(refcoords.cols(), 1) - refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_dx = fe1d_.GradientsReferenceShapeFunctions(Eigen::RowVectorXd::Constant(refcoords.cols(), 1) - refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_dy = fe1d_.GradientsReferenceShapeFunctions(Eigen::RowVectorXd::Constant(refcoords.cols(), 1) - refcoords.row(1));
    for (int i = 0 ; i < refcoords.cols() ; ++i) {
	// Get the gradient of the basis functions associated with the vertices
	result(0, 2*i+0) = sf1d_dx(0, i) * sf1d_y(0, i);
	result(0, 2*i+1) = sf1d_x(0, i) * sf1d_dy(0, i);
	result(1, 2*i+0) = sf1d_dx(1, i) * sf1d_y(0, i);
	result(1, 2*i+1) = sf1d_x(1, i) * sf1d_dy(0, i);
	result(2, 2*i+0) = sf1d_dx(1, i) * sf1d_y(1, i);
	result(2, 2*i+1) = sf1d_x(1, i) * sf1d_dy(1, i);
	result(3, 2*i+0) = sf1d_dx(0, i) * sf1d_y(1, i);
	result(3, 2*i+1) = sf1d_x(0, i) * sf1d_dy(1, i);
	// Get the basis functions associated with the first edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    if (rel_orient_[0] == lf::mesh::Orientation::positive) {
		result(4+j, 2*i+0) = sf1d_dx(2+j, i) * sf1d_y(0, i);
		result(4+j, 2*i+1) = sf1d_x(2+j, i) * sf1d_dy(0, i);
	    }
	    else {
		result(2+degree_-j, 2*i+0) = -sf1df_dx(2+j, i) * sf1d_y(0, i);
		result(2+degree_-j, 2*i+1) = sf1df_x(2+j, i) * sf1d_dy(0, i);
	    }
	}
	// Get the basis functions associated with the second edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    if (rel_orient_[1] == lf::mesh::Orientation::positive) {
		result(3+degree_+j, 2*i+0) = sf1d_dx(1, i) * sf1d_y(2+j, i);
		result(3+degree_+j, 2*i+1) = sf1d_x(1, i) * sf1d_dy(2+j, i);
	    }
	    else {
		result(1+2*degree_-j, 2*i+0) = sf1d_dx(1, i) * sf1df_y(2+j, i);
		result(1+2*degree_-j, 2*i+1) = sf1d_x(1, i) * -sf1df_dy(2+j, i);
	    }
	}
	// Get the basis functions associated with the third edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    if (rel_orient_[2] == lf::mesh::Orientation::positive) {
		result(2+2*degree_+j, 2*i+0) = -sf1df_dx(2+j, i) * sf1d_y(1, i);
		result(2+2*degree_+j, 2*i+1) = sf1df_x(2+j, i) * sf1d_dy(1, i);
	    }
	    else {
		result(3*degree_-j, 2*i+0) = sf1d_dx(2+j, i) * sf1d_y(1, i);
		result(3*degree_-j, 2*i+1) = sf1d_x(2+j, i) * sf1d_dy(1, i);
	    }
	}
	// Get the basis functions associated with the fourth edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    if (rel_orient_[3] == lf::mesh::Orientation::positive) {
		result(1+3*degree_+j, 2*i+0) = sf1d_dx(0, i) * sf1df_y(2+j, i);
		result(1+3*degree_+j, 2*i+1) = sf1d_x(0, i) * -sf1df_dy(2+j, i);
	    }
	    else {
		result(4*degree_-1-j, 2*i+0) = sf1d_dx(0, i) * sf1d_y(2+j, i);
		result(4*degree_-1-j, 2*i+1) = sf1d_x(0, i) * sf1d_dy(2+j, i);
	    }
	}
	// Get the basis functions associated with the interior of the quad
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    for (int k = 0 ; k < degree_-1 ; ++k) {
		result(4*degree_ + (degree_-1)*j +k, 2*i+0) = sf1d_dx(k+2, i) * sf1d_y(j+2, i);
		result(4*degree_ + (degree_-1)*j +k, 2*i+1) = sf1d_x(k+2, i) * sf1d_dy(j+2, i);
	    }
	}
    }
    return result;
  }

  /**
   * @brief Evaluation nodes are the vertices, the Chebyshev nodes of degree p-1
   * on the edges and the corresponding nodes on the quadrilateral
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return eval_nodes_;
  }

  /**
   * @brief (p+1)^2 evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodevals) const override {
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> shape_functions_at_nodes = EvalReferenceShapeFunctions(EvaluationNodes());
    return shape_functions_at_nodes.transpose().fullPivHouseholderQr().solve(nodevals.transpose()).transpose();
  }

 private:
  unsigned degree_;
  Eigen::MatrixXd eval_nodes_;
  FeHPSegment<SCALAR> fe1d_;
  nonstd::span<const lf::mesh::Orientation> rel_orient_;

  Eigen::MatrixXd ComputeEvaluationNodes() const {
    Eigen::MatrixXd nodes(2, (degree_ + 1) * (degree_ + 1));
    const auto cheb = chebyshevNodes(degree_-1);
    // Add the evaluation nodes corresponding to the vertices
    nodes(0, 0) = 0;
    nodes(1, 0) = 0;
    nodes(0, 1) = 1;
    nodes(1, 1) = 0;
    nodes(0, 2) = 1;
    nodes(1, 2) = 1;
    nodes(0, 3) = 0;
    nodes(1, 3) = 1;
    // Add the evaluation nodes corresponding to the edges of the quad
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 4 + i) = cheb[i];
      nodes(1, 4 + i) = 0;
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 3 + degree_ + i) = 1;
      nodes(1, 3 + degree_ + i) = cheb[i];
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 2 + 2*degree_ + i) = 1. - cheb[i];
      nodes(1, 2 + 2*degree_ + i) = 1;
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 1 + 3*degree_ + i) = 0;
      nodes(1, 1 + 3*degree_ + i) = 1. - cheb[i];
    }
    // Add the evaluation nodes corresponding to the interior of the quad
    for (int i = 0; i < degree_ - 1; ++i) {
      for (int j = 0; j < degree_ - 1; ++j) {
        nodes(0, 4 * degree_ + (degree_ - 1) * i + j) = cheb[j];
        nodes(1, 4 * degree_ + (degree_ - 1) * i + j) = cheb[i];
      }
    }
    return nodes;
  }
};

}   // end namespace lf::uscalfe

#endif  // LF_USCALFE_HP_FE_H_

