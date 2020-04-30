#ifndef EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACEHP_H_
#define EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACEHP_H_

#define _USE_MATH_DEFINES

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <cmath>
#include <memory>
#include <vector>


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
 * @brief Lagrangian Finite Elements of arbitrary degreen on segments
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

  FeHPSegment(unsigned degree) : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(), degree_(degree) {
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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override {
    return codim == 0 ? degree_ - 1 : 1;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t,
   * lf::base::sub_idx_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim,
      lf::base::sub_idx_t /*subidx*/) const override {
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
	    return LegendrePoly<SCALAR>::integral(i+2, 2*x-1);
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
	    return 2*LegendrePoly<SCALAR>::eval(i+1, 2*x-1);
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

  Eigen::MatrixXd ComputeEvaluationNodes() const {
      Eigen::MatrixXd nodes(1, degree_+1);
      nodes(0, 0) = 0;
      nodes(0, 1) = 1;
      if (degree_ > 1) {
	  nodes.block(0, 2, 1, degree_-1) = Eigen::VectorXd::LinSpaced(degree_+1, 0, 1).segment(1, degree_-1).transpose();
      }
      return nodes;
  }
};

/**
 * @brief Lagrangian Finite Elements of arbitrary degreen on triangles
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

  FeHPTria(unsigned degree)
      : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(), degree_(degree), eval_nodes_() {
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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override {
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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim,
      lf::base::sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), refcoords.cols());
      // Compute the barycentric coordinate functions
      Eigen::RowVectorXd lambda1 = Eigen::RowVectorXd::Ones(refcoords.cols()) - refcoords.row(0) - refcoords.row(1);
      Eigen::RowVectorXd lambda2 = refcoords.row(0);
      Eigen::RowVectorXd lambda3 = refcoords.row(1);
      // Get the basis functions associated with the vertices
      result.row(0) = lambda1.unaryExpr([&](double x) -> SCALAR { return x; });
      result.row(1) = lambda2.unaryExpr([&](double x) -> SCALAR { return x; });
      result.row(2) = lambda3.unaryExpr([&](double x) -> SCALAR { return x; });
      // Get the basis functions associated with the first edge
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  result.row(i+3) = ((lambda1 + lambda2).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
			     (lambda1.array()/(lambda1+lambda2).array()).unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); })).matrix();
      }
      // Get the basis functions associated with the second edge
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  result.row(i+degree_+2) = ((lambda2 + lambda3).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
			             (lambda2.array()/(lambda2+lambda3).array()).unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); })).matrix();
      }
      // Get the basis functions associated with the third edge
      for (int i = 0 ; i < degree_-1 ; ++i) {
	  result.row(i+2*degree_+1) = ((lambda3 + lambda1).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
			               (lambda3.array()/(lambda3+lambda1).array()).unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); })).matrix();
      }
      // Get the basis functions associated with the interior of the triangle
      if (degree_ > 2) {
	  int idx = 3 * degree_;
	  for (int i = 0 ; i < degree_-2 ; ++i) {
	      for (int j = 0 ; j < degree_-i-2 ; ++j) {
		  result.row(idx) = ((lambda2 + lambda3).unaryExpr([&](double x) -> SCALAR { return std::pow(x, i+2); }).array() *
				     (lambda2.array()/(lambda2+lambda3).array()).unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(i+2, 2*x-1); }) *
				     lambda1.array().unaryExpr([&](double x) -> SCALAR { return LegendrePoly<SCALAR>::integral(j+2, 2*x-1); })).matrix();
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
    Eigen::RowVectorXd lambda1 = Eigen::RowVectorXd::Ones(refcoords.cols()) - refcoords.row(0) - refcoords.row(1);
    Eigen::RowVectorXd lambda2 = refcoords.row(0);
    Eigen::RowVectorXd lambda3 = refcoords.row(1);
    for (int i = 0 ; i < refcoords.cols() ; ++i) {
      // Get the gradient of the basis functions associated with the vertices
      result(0, 2*i+0) = -1;
      result(0, 2*i+1) = -1;
      result(1, 2*i+0) = 1;
      result(1, 2*i+1) = 0;
      result(2, 2*i+0) = 0;
      result(2, 2*i+1) = 1;
      // Get the gradient of the basis functions associated with the first edge
      double lambda1p2 = 1 - refcoords(1, i);
      double lambda1norm = lambda1[i] / lambda1p2;
      for (int j = 0 ; j < degree_-1 ; ++j) {
	  SCALAR leg1inte = LegendrePoly<SCALAR>::integral(j+2, 2*lambda1norm-1);
	  SCALAR leg1eval = LegendrePoly<SCALAR>::eval(j+1, 2*lambda1norm-1);
	  result(j+3, 2*i+0) = -2 * std::pow(lambda1p2, j+1) * leg1eval;
	  result(j+3, 2*i+1) = -(j+2)*std::pow(lambda1p2, j+1)*leg1inte + 2*lambda2[i]*std::pow(lambda1p2, j)*leg1eval;
      }
      // Get the gradient of the basis functions associated with the second edge
      double lambda2p3 = refcoords(0, i) + refcoords(1, i);
      double lambda2norm = lambda2[i] / lambda2p3;
      for (int j = 0 ; j < degree_-1 ; ++j) {
	  SCALAR leg2inte = LegendrePoly<SCALAR>::integral(j+2, 2*lambda2norm-1);
	  SCALAR leg2eval = LegendrePoly<SCALAR>::eval(j+1, 2*lambda2norm-1);
	  result(j+degree_+2, 2*i+0) = (j+2)*std::pow(lambda2p3, j+1)*leg2inte + 2*lambda3[i]*std::pow(lambda2p3, j)*leg2eval;
	  result(j+degree_+2, 2*i+1) = (j+2)*std::pow(lambda2p3, j+1)*leg2inte - 2*lambda2[i]*std::pow(lambda2p3, j)*leg2eval;
      }
      // Get the gradient of the basis functions associated with the third edge
      double lambda3p1 = 1 - refcoords(0, i);
      double lambda3norm = lambda3[i] / lambda3p1;
      for (int j = 0 ; j < degree_-1 ; ++j) {
	  SCALAR leg3inte = LegendrePoly<SCALAR>::integral(j+2, 2*lambda3norm-1);
	  SCALAR leg3eval = LegendrePoly<SCALAR>::eval(j+1, 2*lambda3norm-1);
	  result(j+2*degree_-1, 2*i+0) = -(j+2)*std::pow(lambda3p1, j+1)*leg3inte + 2*lambda3[i]*std::pow(lambda3p1, j)*leg3eval;
	  result(j+2*degree_-1, 2*i+1) = 2 * std::pow(lambda3p1, j+1) * leg3eval;
      }
      // Get the gradient of the basis functions associated with the interior of the triangle
      if (degree_ > 2) {
	  int idx = 3 * degree_;
	  for (int j = 0 ; j < degree_-2 ; ++j) {
	      SCALAR legjinte = LegendrePoly<SCALAR>::integral(j+2, 2*lambda2norm-1);
	      SCALAR legjeval = LegendrePoly<SCALAR>::eval(j+1, 2*lambda2norm-1);
	      for (int k = 0 ; k < degree_-j-2 ; ++k) {
		  SCALAR legkinte = LegendrePoly<SCALAR>::integral(k+2, 2*lambda1[i]-1);
		  SCALAR legkeval = LegendrePoly<SCALAR>::eval(k+1, 2*lambda1[i]-1);
		  result(idx, 2*i+0) = (j+2)*std::pow(lambda2p3, j+1)*legjinte*legkinte + 2*lambda3[i]*std::pow(lambda2p3, j)*legjeval*legkinte - 2*std::pow(lambda2p3, j+2)*legjinte*legkeval;
		  result(idx, 2*i+1) = (j+2)*std::pow(lambda2p3, j+1)*legjinte*legkinte - 2*lambda2[i]*std::pow(lambda2p3, j)*legjeval*legkinte - 2*std::pow(lambda2p3, j+2)*legjinte*legkeval;
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

  Eigen::MatrixXd ComputeEvaluationNodes() const {
    Eigen::MatrixXd eval_nodes(2, (degree_ + 1) * (degree_ + 2) / 2);
    // Add the evaluation nodes corresponding to the vertices of the triangle
    eval_nodes(0, 0) = 0;
    eval_nodes(1, 0) = 0;
    eval_nodes(0, 1) = 1;
    eval_nodes(1, 1) = 0;
    eval_nodes(0, 2) = 0;
    eval_nodes(1, 2) = 1;
    // Add the evaluation nodes corresponding to the edges of the triangle
    for (int i = 0; i < degree_ - 1; ++i) {
      eval_nodes(0, 3 + i) = static_cast<double>(i+1) / (degree_+1);
      eval_nodes(1, 3 + i) = 0;
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      eval_nodes(0, 2 + degree_ + i) = 1. - static_cast<double>(i+1)/(degree_+1);
      eval_nodes(1, 2 + degree_ + i) = static_cast<double>(i+1)/(degree_+1);
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      eval_nodes(0, 1 + 2*degree_ + i) = 0;
      eval_nodes(1, 1 + 2*degree_ + i) = 1. - static_cast<double>(i+1)/(degree_+1);
    }
    // Add the evaluation nodes corresponding to the interior of the triangle
    if (degree_ > 2) {
      int idx = 3 * degree_;
      for (int i = 0; i < degree_ - 2; ++i) {
        for (int j = 0; j < degree_-2-i ; ++j) {
          eval_nodes(0, idx) = static_cast<double>(j+1) / (degree_+1);
          eval_nodes(1, idx) = static_cast<double>(i+1) / (degree_+1);
	  ++idx;
        }
      }
    }
    return eval_nodes;
  }
};

/**
 * @brief Lagrangian Finite Elements of arbitrary degreen on quadrilaterals
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

  FeHPQuad(unsigned degree)
      : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(), degree_(degree), eval_nodes_(), fe1d_(degree) {
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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override {
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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim,
      lf::base::sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(NumRefShapeFunctions(), refcoords.cols());
    // Compute the 1D shape functions at the x and y coordinates
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_x = fe1d_.EvalReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_y = fe1d_.EvalReferenceShapeFunctions(refcoords.row(1));
    // Get the basis functions associated with the vertices
    result.row(0) = (sf1d_x.row(0).array() * sf1d_y.row(0).array()).matrix();
    result.row(1) = (sf1d_x.row(1).array() * sf1d_y.row(0).array()).matrix();
    result.row(2) = (sf1d_x.row(1).array() * sf1d_y.row(1).array()).matrix();
    result.row(3) = (sf1d_x.row(0).array() * sf1d_y.row(1).array()).matrix();
    // Get the basis functions associated with the first edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	result.row(4+i) = (sf1d_x.row(i+2).array() * sf1d_y.row(0).array()).matrix();
    }
    // Get the basis functions associated with the second edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	result.row(3+degree_+i) = (sf1d_x.row(1).array() * sf1d_y.row(i+2).array()).matrix();
    }
    // Get the basis functions associated with the third edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	//result.row(2+2*degree_+i) = (sf1d_x.row(degree_-i).array() * sf1d_y.row(1).array()).matrix();
	result.row(2+2*degree_+i) = (sf1d_x.row(i+2).array() * sf1d_y.row(1).array()).matrix();
    }
    // Get the basis functions associated with the fourth edge
    for (int i = 0 ; i < degree_-1 ; ++i) {
	//result.row(1+3*degree_+i) = (sf1d_x.row(0).array() * sf1d_y.row(degree_-i).array()).matrix();
	result.row(1+3*degree_+i) = (sf1d_x.row(0).array() * sf1d_y.row(i+2).array()).matrix();
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
	    result(4+j, 2*i+0) = sf1d_dx(j+2, i) * sf1d_y(0, i);
	    result(4+j, 2*i+1) = sf1d_x(j+2, i) * sf1d_dy(0, i);
	}
	// Get the basis functions associated with the second edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    result(3+degree_+j, 2*i+0) = sf1d_dx(1, i) * sf1d_y(j+2, i);
	    result(3+degree_+j, 2*i+1) = sf1d_x(1, i) * sf1d_dy(j+2, i);
	}
	// Get the basis functions associated with the third edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    //result(2+2*degree_+j, 2*i+0) = sf1d_dx(degree_-j, i) * sf1d_y(1, i);
	    //result(2+2*degree_+j, 2*i+1) = sf1d_x(degree_-j, i) * sf1d_dy(1, i);
	    result(2+2*degree_+j, 2*i+0) = sf1d_dx(j+2, i) * sf1d_y(1, i);
	    result(2+2*degree_+j, 2*i+1) = sf1d_x(j+2, i) * sf1d_dy(1, i);
	}
	// Get the basis functions associated with the fourth edge
	for (int j = 0 ; j < degree_-1 ; ++j) {
	    //result(1+3*degree_+j, 2*i+0) = sf1d_dx(0, i) * sf1d_y(degree_-j, i);
	    //result(1+3*degree_+j, 2*i+1) = sf1d_x(0, i) * sf1d_dy(degree_-j, i);
	    result(1+3*degree_+j, 2*i+0) = sf1d_dx(0, i) * sf1d_y(j+2, i);
	    result(1+3*degree_+j, 2*i+1) = sf1d_x(0, i) * sf1d_dy(j+2, i);
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

  Eigen::MatrixXd ComputeEvaluationNodes() const {
    Eigen::MatrixXd nodes(2, (degree_ + 1) * (degree_ + 1));
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
      nodes(0, 4 + i) = static_cast<double>(i+1) / (degree_+1);
      nodes(1, 4 + i) = 0;
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 3 + degree_ + i) = 1;
      nodes(1, 3 + degree_ + i) = static_cast<double>(i+1) / (degree_+1);
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 2 + 2*degree_ + i) = 1. - static_cast<double>(i+1)/(degree_+1);
      nodes(1, 2 + 2*degree_ + i) = 1;
    }
    for (int i = 0; i < degree_ - 1; ++i) {
      nodes(0, 1 + 3*degree_ + i) = 0;
      nodes(1, 1 + 3*degree_ + i) = 1. - static_cast<double>(i+1)/(degree_+1);
    }
    // Add the evaluation nodes corresponding to the interior of the quad
    for (int i = 0; i < degree_ - 1; ++i) {
      for (int j = 0; j < degree_ - 1; ++j) {
        nodes(0, 4 * degree_ + (degree_ - 1) * i + j) = static_cast<double>(j+1) / (degree_+1);
        nodes(1, 4 * degree_ + (degree_ - 1) * i + j) = static_cast<double>(i+1) / (degree_+1);
      }
    }
    return nodes;
  }
};

/**
 * @brief Lagrangian Finite Element Space of arbitrary degree
 */
template <typename SCALAR>
class FeSpaceHP : public lf::uscalfe::UniformScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  FeSpaceHP() = delete;
  FeSpaceHP(const FeSpaceHP &) = delete;
  FeSpaceHP(FeSpaceHP &&) noexcept = default;
  FeSpaceHP &operator=(const FeSpaceHP &) = delete;
  FeSpaceHP &operator=(FeSpaceHP &&) noexcept = default;

  /**
   * @brief Constructor: Sets up the dof handler
   * @param mesh_p A shared pointer to the underlying mesh (immutable)
   * @param N The polynomial degree of the Finite Element Space
   */
  explicit FeSpaceHP(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, unsigned N)
      : lf::uscalfe::UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<FeHPTria<SCALAR>>(N),
            std::make_shared<FeHPQuad<SCALAR>>(N),
            std::make_shared<FeHPSegment<SCALAR>>(N),
            std::make_shared<lf::uscalfe::FeLagrangePoint<SCALAR>>(N)) {
	// This constructor does very unorthodox things, but there is no easy way around it
	// with the current state of LehrFEM++
	// Permute the edge dofs such that they have a global instead of a local ordering
	const auto& dofh = lf::uscalfe::UniformScalarFESpace<SCALAR>::LocGlobMap();
	for (auto cell : mesh_p->Entities(0)) {
	    const auto dofidxs = dofh.GlobalDofIndices(*cell);
	    // Extreamly evil but necessary as UniformFEDofHandler does not expose its dofs_ array
	    lf::assemble::gdof_idx_t* edge_dofidx = const_cast<lf::assemble::gdof_idx_t*>(dofidxs.data());
	    edge_dofidx += cell->RefEl().NumSubEntities(2);
	    const auto orient = cell->RelativeOrientations();
	    const auto edges = cell->SubEntities(1);
	    // Iterate over all edges and reverse their dofs if orient[i] is negative
	    for (int i = 0 ; i < edges.size() ; ++i) {
		const auto num_edge_dofs = dofh.NumInteriorDofs(*(edges[i]));
		if (orient[i] == lf::mesh::Orientation::negative) {
		    std::reverse(edge_dofidx, edge_dofidx+num_edge_dofs);
		}
		edge_dofidx += num_edge_dofs;
	    }
	}
    }

  ~FeSpaceHP() override = default;
};

#endif  // EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACEHP_H_
