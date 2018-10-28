#ifndef LF_LOCCOMPELLBVP
#define LF_LOCCOMPELLBVP
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Classes taking care of local computations for scalar 2nd-order
 * elliptic BVPs
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <lf/quad/quad.h>
#include <iostream>
#include "lagr_fe.h"

namespace lf::fe {
/** @brief Computing the element matrix for the (negative) Laplacian
 *         and linear finite elements.
 *
 * The main purpose of this class is to compute the element matrix for
 * the Laplacian on affine triangles or bilinearly mapped quadrilaterals.
 * These element matrices are provided by the `Eval()` method.
 *
 * @note the `Eval()` method will always return a _reference_ to a 4x4 matrix
 * also for triangles. In this case the last row and column must be ignored.
 *
 * This class complies with the requirements for the type `ELEM_MAT_COMP`
 * given as a template parameter to define an incarnation of the function
 * AssembleMatrixLocally().
 */
class LinearFELaplaceElementMatrix {
 public:
  using elem_mat_t = Eigen::Matrix<double, 4, 4>;
  using ElemMat = const elem_mat_t;

  /**
   * @brief Idle constructor
   */
  LinearFELaplaceElementMatrix() {}
  /**
   * @brief All cells are considered active in the default implementation
   */
  virtual bool isActive(const lf::mesh::Entity &cell) { return true; }
  /*
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a 4x4 matrix, containing the element matrix. The bottom row/column
   *         is not used in the case of a triangle.
   */
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  /// quadrature points on reference quadrilateral
  const double kSqrt3 = 1.0 / std::sqrt(3.0);
  const double kZeta_0 = 0.5 - 0.5 * kSqrt3;
  const double kZeta_1 = 0.5 + 0.5 * kSqrt3;
  const std::array<Eigen::Vector2d, 4> kQuadPoints{
      Eigen::Vector2d{kZeta_0, kZeta_0}, Eigen::Vector2d{kZeta_0, kZeta_1},
      Eigen::Vector2d{kZeta_1, kZeta_0}, Eigen::Vector2d{kZeta_1, kZeta_1}};
  // Gradients of refrence shape functions in rows of a matrix
  static Eigen::Matrix<double, 4, 2> DervRefShapFncts(
      const Eigen::Vector2d &xh);
  // Barycenter of triangle
  const Eigen::Matrix<double, 2, 1> kTriabc{
      (Eigen::Matrix<double, 2, 1>() << (1.0 / 3), (1.0 / 3)).finished()};
  // Point in reference square
  const Eigen::Matrix<double, 2, 1> kQuadpt{
      (Eigen::Matrix<double, 2, 1>() << 0.7, 0.3).finished()};

 public:
  /*
   * @brief static variable for controling debugging output
   */
  static unsigned int dbg_ctrl;
  static const unsigned int dbg_det = 1;
  static const unsigned int dbg_locmat = 2;
  static const unsigned int dbg_J = 4;
  static const unsigned int dbg_geo = 8;
};

/** @brief Class for computation of local load vector for linear finite
 * elements.
 *
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<double(const Eigen::Vector2d &)>, which supplies
 *         the source function
 *
 * Computation is based on vertex based quadrature
 *
 * @note The element vector returned by the `Eval()` method will always
 * have length 4 also for triangles.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 */
template <typename SCALAR, typename FUNCTOR>
class LinearFELocalLoadVector {
 public:
  using elem_vec_t = Eigen::Matrix<SCALAR, 4, 1>;
  using ElemVec = const elem_vec_t;

  /** @brief Constructor storing the right hand side function */
  LinearFELocalLoadVector(FUNCTOR f) : f_(f) {}
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity &cell) { return true; }
  /*
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   *
   * The implementation uses simple vertex based quadrature and an approximation
   * of the volume of a cell just using the integration element at the
   * barycenter.
   */
  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  FUNCTOR f_;

 public:
  /*
   * @brief static variable for controling debugging output
   */
  static unsigned int dbg_ctrl;
  static const unsigned int dbg_locvec = 1;
  static const unsigned int dbg_geo = 2;
};

template <typename SCALAR, typename FUNCTOR>
unsigned int LinearFELocalLoadVector<SCALAR, FUNCTOR>::dbg_ctrl = 0;

template <typename SCALAR, typename FUNCTOR>
typename LinearFELocalLoadVector<SCALAR, FUNCTOR>::ElemVec
LinearFELocalLoadVector<SCALAR, FUNCTOR>::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  const Eigen::MatrixXd &ref_el_corners(ref_el.NodeCoords());
  const Eigen::MatrixXd vertices{geo_ptr->Global(ref_el_corners)};
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_geo,
                    std::cout << ref_el << ", shape = \n"
                              << vertices << std::endl);
  const double area = lf::geometry::Volume(*geo_ptr);

  // Vector for returning element vector

  elem_vec_t elem_vec = elem_vec_t::Zero();
  // get function values in the vertices
  for (int k = 0; k < num_nodes; k++) {
    elem_vec[k] = area * f_(vertices.col(k)) / num_nodes;
  }
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_locvec,
      std::cout << "element vector = " << elem_vec.head(num_nodes).transpose()
                << std::endl);
  return elem_vec;
}

/**
 * @brief Class for local quadrature based computations for Lagrangian finite
 * elements
 *
 *
 */
template <typename DIFF_COEFF, typename REACTION_COEFF>
class LagrangeFEEllBVPElementMatrix {
 public:
  /**
   * @brief type of returned element matrix
   */
  using elem_mat_t = Eigen::MatrixXd;
  using ElemMat = const elem_mat_t;
  /*
   * @brief Constructor: cell-independent precomputations
   *
   * @param fe_trie finite element to be used on triangles
   * @param fe_quad finite element for quadrilaterals
   *
   * The two parametric finite elements for triangles and quadrilaterals have to
   * be compatible in the sense that they all assign exactly one reference shape
   * function to each vertex and the same number of interior shape functions
   * to each edge
   */
  LagrangeFEEllBVPElementMatrix(
      const ScalarReferenceFiniteElement<double> &fe_tria,
      const ScalarReferenceFiniteElement<double> &fe_quad, DIFF_COEFF alpha,
      REACTION_COEFF gamma);
  /**
   * @brief All cells are considered active in the default implementation
   */
  virtual bool isActive(const lf::mesh::Entity &cell) { return true; }
  /*
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a small dense, containing the element matrix.
   */
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  /**
   * @brief functors providing coefficient functions
   */
  DIFF_COEFF alpha_;
  REACTION_COEFF gamma_;

  /**
   * @brief references to finite elements
   */
  const ScalarReferenceFiniteElement<double> &fe_tria_, &fe_quad_;
  /**
   * @brief Numbers of local shape functions
   */
  size_type Nrsf_tria_, Nrsf_quad_;
  /**
   * @brief Numbers of quadrature points
   */
  size_type Nqp_tria_, Nqp_quad_;
  /**
   * @brief quadrature rules for both admissible types of cells
   */
  lf::quad::QuadRule qr_tria_, qr_quad_;
  /**
   * @brief Matrix of values of all reference shape functions at all quadrature
   * points
   *
   * The rows correspond to the RSFs, the columns to the quadrature points
   */
  Eigen::MatrixXd rsf_quadpoints_tria_, rsf_quadpoints_quad_;
  /**
   * @brief array of matrices whose columns contain the gradients of all RFSs in
   *        a single quadrature point
   */
  std::vector<Eigen::Matrix<double, 2, Eigen::Dynamic>> grad_quadpoint_tria_,
      grad_quadpoint_quad_;

 public:
  /** @brief output control variable
   */
  static unsigned int ctrl_;
  static const unsigned int kout_qr = 1;
  static const unsigned int kout_rsfvals = 2;
  static const unsigned int kout_gradvals = 4;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_locmat = 16;
};

template <typename DIFF_COEFF, typename REACTION_COEFF>
unsigned int LagrangeFEEllBVPElementMatrix<DIFF_COEFF, REACTION_COEFF>::ctrl_ =
    0;

template <typename DIFF_COEFF, typename REACTION_COEFF>
LagrangeFEEllBVPElementMatrix<DIFF_COEFF, REACTION_COEFF>::
    LagrangeFEEllBVPElementMatrix(
        const ScalarReferenceFiniteElement<double> &fe_tria,
        const ScalarReferenceFiniteElement<double> &fe_quad, DIFF_COEFF alpha,
        REACTION_COEFF gamma)
    : fe_tria_(fe_tria), fe_quad_(fe_quad), alpha_(alpha), gamma_(gamma) {
  LF_ASSERT_MSG((fe_tria_.Dimension() == 2) && (fe_quad_.Dimension() == 2),
                "Implemented only in 2D!");
  LF_ASSERT_MSG((fe_tria_.RefEl() == lf::base::RefEl::kTria()) ||
                    (fe_quad_.RefEl() == lf::base::RefEl::kQuad()),
                "Unexpected type of reference cell");
  LF_ASSERT_MSG((fe_tria_.NumRefShapeFunctions(2) == 1) &&
                    (fe_quad_.NumRefShapeFunctions(2) == 1),
                "Exactly one shape function must be assigned to each vertex");
  LF_ASSERT_MSG(
      (fe_tria_.NumRefShapeFunctions(1) == fe_quad_.NumRefShapeFunctions(1)),
      "#RSF mismatch on edges " << fe_tria_.NumRefShapeFunctions(1) << " <-> "
                                << fe_quad_.NumRefShapeFunctions(1));

  // Maximal order of both finite elements
  const unsigned int poly_order = std::max(fe_tria_.order(), fe_quad_.order());

  // Obtain quadrature rules for both triangles and quadrilaterals
  // We choose the order twice the polynomial degree of the finite element
  // space. This is slightly more than required for an admissible variational
  // crime.
  lf::quad::quadOrder_t quad_order = 2 * poly_order;

  {
    // Preprocessing for triangles
    Nrsf_tria_ = fe_tria_.NumRefShapeFunctions();

    // Fetch suitable predefined quadrature rule
    qr_tria_ = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), quad_order);
    Nqp_tria_ = qr_tria_.NumPoints();
    SWITCHEDSTATEMENT(ctrl_, kout_qr,
                      std::cout << "LagrEM(Tria): " << qr_tria_ << std::endl);

    // Obtain value of reference shape functions in all quadrature points
    const std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> rsf_val{
        fe_tria_.EvalReferenceShapeFunctions(qr_tria_.Points())};
    LF_ASSERT_MSG(Nrsf_tria_ == rsf_val.size(),
                  "Mismatch in length of value vector " << Nrsf_tria_ << " <-> "
                                                        << rsf_val.size());
    // Store the values for reuse for all cells
    rsf_quadpoints_tria_.resize(Nrsf_tria_, Nqp_tria_);

    for (int i = 0; i < Nrsf_tria_; ++i) {
      LF_ASSERT_MSG(
          (rsf_val[i].cols() == Nqp_tria_),
          "Length mismatch " << rsf_val[i].cols() << " <-> " << Nqp_tria_);
      for (int j = 0; j < Nqp_tria_; ++j) {
        rsf_quadpoints_tria_(i, j) = rsf_val[i][j];
      }
    }
    SWITCHEDSTATEMENT(ctrl_, kout_rsfvals,
                      std::cout << "LagrEM(Tria): values of RSFs\n"
                                << rsf_quadpoints_tria_ << std::endl);

    // Store the gradients for the reference shape functions
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> rsf_grad{
        fe_tria_.GradientsReferenceShapeFunctions(qr_tria_.Points())};
    LF_ASSERT_MSG(Nrsf_tria_ == rsf_grad.size(),
                  "Mismatch in length of gradient vector"
                      << Nrsf_tria_ << " <-> " << rsf_grad.size());

    // Store the gradients internally in an array of stacked vectors
    grad_quadpoint_tria_.resize(
        Nqp_tria_, Eigen::Matrix<double, 2, Eigen::Dynamic>(2, Nrsf_tria_));
    for (int i = 0; i < Nqp_tria_; ++i) {
      for (int j = 0; j < Nrsf_tria_; ++j) {
        grad_quadpoint_tria_[i].col(j) = rsf_grad[j].col(i);
      }
    }
    SWITCHEDSTATEMENT(ctrl_, kout_gradvals,
                      std::cout << "LagrEM(Tria): gradients:" << std::endl;
                      for (int i = 0; i < Nqp_tria_; ++i) {
                        std::cout << "QP " << i << " = \n"
                                  << grad_quadpoint_tria_[i] << std::endl;
                      });
  }

  {
    // Preprocessing for quadrilaterals
    Nrsf_quad_ = fe_quad_.NumRefShapeFunctions();

    // Fetch suitable predefined quadrature rule
    qr_quad_ = lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), quad_order);
    Nqp_quad_ = qr_quad_.NumPoints();
    SWITCHEDSTATEMENT(ctrl_, kout_qr,
                      std::cout << "LagrEM(Quad): " << qr_quad_ << std::endl);

    // Obtain value of reference shape functions in all quadrature points
    const std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> rsf_val{
        fe_quad_.EvalReferenceShapeFunctions(qr_quad_.Points())};
    LF_ASSERT_MSG(Nrsf_quad_ == rsf_val.size(),
                  "Mismatch in length of value vector " << Nrsf_quad_ << " <-> "
                                                        << rsf_val.size());
    // Store the values for reuse for all cells
    rsf_quadpoints_quad_.resize(Nrsf_quad_, Nqp_quad_);

    for (int i = 0; i < Nrsf_quad_; ++i) {
      LF_ASSERT_MSG(
          (rsf_val[i].cols() == Nqp_quad_),
          "Length mismatch " << rsf_val[i].cols() << " <-> " << Nqp_quad_);
      for (int j = 0; j < Nqp_quad_; ++j) {
        rsf_quadpoints_quad_(i, j) = rsf_val[i][j];
      }
    }
    SWITCHEDSTATEMENT(ctrl_, kout_rsfvals,
                      std::cout << "LagrEM(Quad): values of RSFs\n"
                                << rsf_quadpoints_quad_ << std::endl);

    // Store the gradients for the reference shape functions
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> rsf_grad{
        fe_quad_.GradientsReferenceShapeFunctions(qr_quad_.Points())};
    LF_ASSERT_MSG(Nrsf_quad_ == rsf_grad.size(),
                  "Mismatch in length of gradient vector"
                      << Nrsf_quad_ << " <-> " << rsf_grad.size());

    // Store the gradients internally in an array of stacked vectors
    grad_quadpoint_quad_.resize(
        Nqp_quad_, Eigen::Matrix<double, 2, Eigen::Dynamic>(2, Nrsf_quad_));
    for (int i = 0; i < Nqp_quad_; ++i) {
      for (int j = 0; j < Nrsf_quad_; ++j) {
        grad_quadpoint_quad_[i].col(j) = rsf_grad[j].col(i);
      }
    }
    SWITCHEDSTATEMENT(ctrl_, kout_gradvals,
                      std::cout << "LagrEM(Quad): gradients:" << std::endl;
                      for (int i = 0; i < Nqp_quad_; ++i) {
                        std::cout << "QP " << i << " = \n"
                                  << grad_quadpoint_quad_[i] << std::endl;
                      });
  }
}  // end constructor LagrangeFEEllBVPElementMatrix<DIFF_COEFF, REACTION_COEFF>

template <typename DIFF_COEFF, typename REACTION_COEFF>
typename LagrangeFEEllBVPElementMatrix<DIFF_COEFF, REACTION_COEFF>::ElemMat
LagrangeFEEllBVPElementMatrix<DIFF_COEFF, REACTION_COEFF>::Eval(
    const lf::mesh::Entity &cell) {
  // Type for diffusion coefficient
  using diff_coeff_t = decltype(alpha_(Eigen::Vector2d::Zero()));
  // Type for reaction coefficient
  using reac_coeff_t = decltype(gamma_(Eigen::Vector2d::Zero()));

  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SWITCHEDSTATEMENT(ctrl_, kout_cell,
                    std::cout << ref_el << ", shape = \n"
                              << geo_ptr->Global(ref_el.NodeCoords())
                              << std::endl);

  // Computations differ depending on the type of the cell
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      // Quadrature points
      const Eigen::Matrix quadpts(qr_tria_.Points());
      // Quadrature weights
      const Eigen::VectorXd quadw(qr_tria_.Weights());
      // Obtain the metric factors for the quadrature points
      const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(quadpts));
      // Fetch the transformation matrices for the gradients
      const Eigen::MatrixXd JinvT(geo_ptr->JacobianInverseGramian(quadpts));
      // Compute the coefficients at quadrature points
      std::vector<diff_coeff_t> alphavals{};
      std::vector<reac_coeff_t> gammavals{};
      for (int i = 0; i < Nqp_tria_; ++i) {
        alphavals.push_back(alpha_(quadpts.col(i)));
        gammavals.push_back(gamma_(quadpts.col(i)));
      }

      // Element matrix
      Eigen::MatrixXd mat(Nrsf_tria_, Nrsf_tria_);
      mat.setZero();

      // Loop over quadrature points
      for (int k = 0; k < Nqp_tria_; ++k) {
        const double w = quadw[k] * determinants[k];
        // Transformed gradients
        const auto trf_grad(JinvT.block(0, 2 * k, 2, 2) *
                            grad_quadpoint_tria_[k]);
        // Transformed gradients multiplied with coefficient
        const auto alpha_trf_grad(alphavals[k] * trf_grad);
        mat += w * (alpha_trf_grad.transpose() * trf_grad +
                    (gammavals[k] * rsf_quadpoints_tria_.col(k)) *
                        (rsf_quadpoints_tria_.col(k).transpose()));
      }
      return mat;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      // Quadrature points
      const Eigen::Matrix quadpts(qr_quad_.Points());
      // Quadrature weights
      const Eigen::VectorXd quadw(qr_quad_.Weights());
      // Obtain the metric factors for the quadrature points
      const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(quadpts));
      // Fetch the transformation matrices for the gradients
      const Eigen::MatrixXd JinvT(geo_ptr->JacobianInverseGramian(quadpts));
      // Compute the coefficients at quadrature points
      std::vector<diff_coeff_t> alphavals{};
      std::vector<reac_coeff_t> gammavals{};
      for (int i = 0; i < Nqp_quad_; ++i) {
        alphavals.push_back(alpha_(quadpts.col(i)));
        gammavals.push_back(gamma_(quadpts.col(i)));
      }

      // Element matrix
      Eigen::MatrixXd mat(Nrsf_quad_, Nrsf_quad_);
      mat.setZero();

      // Loop over quadrature points
      for (int k = 0; k < Nqp_quad_; ++k) {
        const double w = quadw[k] * determinants[k];
        // Transformed gradients
        const auto trf_grad(JinvT.block(0, 2 * k, 2, 2) *
                            grad_quadpoint_quad_[k]);
        // Transformed gradients multiplied with coefficient
        const auto alpha_trf_grad(alphavals[k] * trf_grad);
        mat += w * (alpha_trf_grad.transpose() * trf_grad +
                    (gammavals[k] * rsf_quadpoints_quad_.col(k)) *
                        (rsf_quadpoints_quad_.col(k).transpose()));
      }
      return mat;
      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch
  return Eigen::MatrixXd(0,0);
}
}  // namespace lf::fe

#endif
