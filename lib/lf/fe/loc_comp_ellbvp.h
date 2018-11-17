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
  LinearFELaplaceElementMatrix() = default;
  /**
   * @brief All cells are considered active in the default implementation
   */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
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
  explicit LinearFELocalLoadVector(FUNCTOR f) : f_(f) {}
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
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
  // World coordinates of vertices
  // const Eigen::MatrixXd vertices{geo_ptr->Global(ref_el_corners)};
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_geo,
                    std::cout << ref_el << ", shape = \n"
                              << (geo_ptr->Global(ref_el_corners))
                              << std::endl);
  // Midpoints of edges in the reference cell
  Eigen::MatrixXd ref_mp(2, 4);
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      // clang-format off
      ref_mp << 0.5,0.5,0.0,-1.0,
                0.0,0.5,0.5,-1.0;
      // clang-format on
      break;
    }
    case lf::base::RefEl::kQuad(): {
      // clang-format off
      ref_mp << 0.5,1,0.5,0.0,
                0.0,0.5,1.0,0.5;
      // clang-format on
      break;
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal entity type!");
      break;
    }
  }  // end switch
  // Midpoints of edges in world coordinates
  const Eigen::MatrixXd mp(geo_ptr->Global(ref_mp));

  const double area = lf::geometry::Volume(*geo_ptr);

  // Vector for returning element vector
  elem_vec_t elem_vec = elem_vec_t::Zero();
  // Run over the midpoints of edges and fetch values of the source function
  // there
  for (int k = 0; k < num_nodes; k++) {
    const auto fval_half = 0.5 * f_(mp.col(k));
    elem_vec[k] += fval_half;
    elem_vec[(k + 1) % num_nodes] += fval_half;
  }
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_locvec,
      std::cout << "element vector = " << elem_vec.head(num_nodes).transpose()
                << std::endl);
  return (area / num_nodes) * elem_vec;
}

/**
 * @brief Preparing reusable information for parametric Lagrangian finite
 * elements
 *
 * For _uniform_ _parametric_ Lagrangian finite element spaces local
 * computations applied to all celss of a mesh often rely on a single quadrature
 * rule for a particular cell type. Therefore the values and gradients of
 * reference shape functions can be precomputed at reference quadrature points
 * and stored for reuse on every cell. This is the purpose of this class.
 */
class LocCompLagrFEPreprocessor {
 protected:
  /* @brief Constructor: initializing local data
   *
   * @param fe_trie finite element to be used on triangles
   * @param fe_quad finite element for quadrilaterals
   * @param loc_quad_order desired order of local quadrature, default value = 0.
   *        If = 0, the quadrature order is determined from the polynomial
   *        degree of the reference shape functions.
   *
   * The two parametric finite elements for triangles and quadrilaterals have to
   * be compatible in the sense that they all assign exactly one reference shape
   * function to each vertex and the same number of interior shape functions
   * to each edge.
   *
   * The constructor sets up local quadrature rules and fetches the values and
   * gradients of references local shape functions in the quadrature points (on
   * the reference element)
   */
  LocCompLagrFEPreprocessor(const ScalarReferenceFiniteElement<double> &fe_tria,
                            const ScalarReferenceFiniteElement<double> &fe_quad,
                                lf::quad::quadOrder_t loc_quad_order = 0);

  /** @brief type-dependent query of quadrature points */
  Eigen::MatrixXd qr_points(lf::base::RefEl ref_el) {
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        return qr_tria_.Points();
      }
      case lf::base::RefEl::kQuad(): {
        return qr_quad_.Points();
      }
      default: { LF_VERIFY_MSG(false, "Illegal type " << ref_el); }
    }
    return Eigen::MatrixXd(0, 0);
  }
  /** @brief type-dependent query of quadrature weights */
  Eigen::VectorXd qr_weights(lf::base::RefEl ref_el) {
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        return qr_tria_.Weights();
      }
      case lf::base::RefEl::kQuad(): {
        return qr_quad_.Weights();
      }
      default: { LF_VERIFY_MSG(false, "Illegal type " << ref_el); }
    }
    return Eigen::VectorXd(0);
  }
  /** @brief type-dependent query of number of quadrature points */
  size_type qr_num_pts(lf::base::RefEl ref_el) {
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        return Nqp_tria_;
      }
      case lf::base::RefEl::kQuad(): {
        return Nqp_quad_;
      }
      default: { LF_VERIFY_MSG(false, "Illegal type " << ref_el); }
    }
    return 0;
  }
  /** @brief type-dependent total number of local shape functions */
  size_type num_rsf(lf::base::RefEl ref_el) {
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        return Nrsf_tria_;
      }
      case lf::base::RefEl::kQuad(): {
        return Nrsf_quad_;
      }
      default: { LF_VERIFY_MSG(false, "Illegal type " << ref_el); }
    }
    return 0;
  }
  /** @brief Matrix of values of reference shape functions at theorem
   *        quadrature points in the reference element
   * @param type of reference element (triangle or quadrilateral)
   * @return an N_rsf x N_qpt matrix, each row corresponding to a reference
   *         shape function, each column to a quadrature point.
   */
  const Eigen::MatrixXd &rsf_at_quadpts(lf::base::RefEl ref_el) {
    if (ref_el == lf::base::RefEl::kQuad()) {
      return rsf_quadpoints_quad_;
    } else
      return rsf_quadpoints_tria_;
  }
  /** @brief Vector of packed gradients of reference shape funtions at
   *         reference quadrature points
   * @param type of reference element (triangle or quadrilateral)
   * @return array of matrices whose columns contain the gradients of all RFSs
   * in a single quadrature point. The length of the array agrees with the
   * number of quadrature points
   */
  const std::vector<Eigen::Matrix<double, 2, Eigen::Dynamic>>
      &grad_rsf_at_quadpts(lf::base::RefEl ref_el) {
    if (ref_el == lf::base::RefEl::kQuad()) {
      return grad_quadpoint_quad_;
    } else
      return grad_quadpoint_tria_;
  }

  /**
   * @brief references to arrangement of local shape functions
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
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_qr = 1;
  static const unsigned int kout_rsfvals = 2;
  static const unsigned int kout_gradvals = 4;
};

/**
 * @brief Class for local quadrature based computations for Lagrangian finite
 * elements
 *
 * @tparam SCALAR type for the entries of the element matrices
 * @tparam DIFF_COEFF a functor providing point evaluation for the diffusion
 * tensor
 * @tparam REACTION_COEFF a functor for point evaluation of the reaction
 * coefficient
 *
 * ### Template parameter type requirements
 *
 * - SCALAR must be a field type like `std::complex<double>`
 * - DIFF_COEFF must comply with `std::function<T(Eigen::Vector2d)>`, where `T`
 * is either a SCALAR compatible type of a matrix type like
 * `Eigen::Matrix<SCALAR,...>`.
 * - REACTION_COEFF must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 *
 * This class complies with the type requirements for the template argument
 * ELEM_MAT_COMP of the function lf::assemble::AssembleMatrixLocally().
 */
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
class LagrangeFEEllBVPElementMatrix : public LocCompLagrFEPreprocessor {
 public:
  /**
   * @brief type of returned element matrix
   */
  using elem_mat_t = Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>;
  using ElemMat = const elem_mat_t;
  /** Type for diffusion coefficient */
  using diff_coeff_t =
      typename std::invoke_result<DIFF_COEFF, Eigen::Vector2d>::type;
  /** Type for reaction coefficient */
  using reac_coeff_t =
      typename std::invoke_result<REACTION_COEFF, Eigen::Vector2d>::type;

  /*
   * @brief Constructor: cell-independent precomputations
   *
   * @param fe_trie finite element to be used on triangles
   * @param fe_quad finite element for quadrilaterals
   *
   * The two parametric finite elements for triangles and quadrilaterals have to
   * be compatible in the sense that they all assign exactly one reference shape
   * function to each vertex and the same number of interior shape functions
   * to each edge.
   *
   * For the sake of efficiency, this constructor precomputes the values of
   * reference shape functions and their gradients at the quadrature points in
   * the reference element.
   */
  LagrangeFEEllBVPElementMatrix(
      const ScalarReferenceFiniteElement<double> &fe_tria,
      const ScalarReferenceFiniteElement<double> &fe_quad, DIFF_COEFF alpha,
      REACTION_COEFF gamma);
  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /*
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a small dense, containing the element matrix.
   *
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials o degree 2p.
   */
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  /**
   * @brief functors providing coefficient functions
   */
  /** @{ */
  DIFF_COEFF alpha_;
  REACTION_COEFF gamma_;
  /** @} */
 public:
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_locmat = 16;
};

template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
unsigned int
    LagrangeFEEllBVPElementMatrix<SCALAR, DIFF_COEFF, REACTION_COEFF>::ctrl_ =
        0;

template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
LagrangeFEEllBVPElementMatrix<SCALAR, DIFF_COEFF, REACTION_COEFF>::
    LagrangeFEEllBVPElementMatrix(
        const ScalarReferenceFiniteElement<double> &fe_tria,
        const ScalarReferenceFiniteElement<double> &fe_quad, DIFF_COEFF alpha,
        REACTION_COEFF gamma)
    : LocCompLagrFEPreprocessor(fe_tria, fe_quad),
      alpha_(alpha),
      gamma_(gamma) {}

template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
typename LagrangeFEEllBVPElementMatrix<SCALAR, DIFF_COEFF,
                                       REACTION_COEFF>::ElemMat
LagrangeFEEllBVPElementMatrix<SCALAR, DIFF_COEFF, REACTION_COEFF>::Eval(
    const lf::mesh::Entity &cell) {
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
      // Quadrature points in actual cell
      const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr_tria_.Points()));
      // Obtain the metric factors for the quadrature points
      const Eigen::VectorXd determinants(
          geo_ptr->IntegrationElement(qr_tria_.Points()));
      // Fetch the transformation matrices for the gradients
      const Eigen::MatrixXd JinvT(
          geo_ptr->JacobianInverseGramian(qr_tria_.Points()));

      // Element matrix
      elem_mat_t mat(Nrsf_tria_, Nrsf_tria_);
      mat.setZero();

      // Loop over quadrature points
      for (int k = 0; k < Nqp_tria_; ++k) {
        // Evaluate diffusion and reaction coefficient
        // at quadrature points in actual cell
        const auto alphaval = alpha_(mapped_qpts.col(k));
        const auto gammaval = gamma_(mapped_qpts.col(k));

        const double w = qr_tria_.Weights()[k] * determinants[k];
        // Transformed gradients
        const auto trf_grad(JinvT.block(0, 2 * k, 2, 2) *
                            grad_quadpoint_tria_[k]);
        // Transformed gradients multiplied with coefficient
        const auto alpha_trf_grad(alphaval * trf_grad);
        mat += w * (alpha_trf_grad.transpose() * trf_grad +
                    (gammaval * rsf_quadpoints_tria_.col(k)) *
                        (rsf_quadpoints_tria_.col(k).transpose()));
      }
      return mat;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      // Quadrature points in actual cell
      const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr_quad_.Points()));
      // Obtain the metric factors for the quadrature points
      const Eigen::VectorXd determinants(
          geo_ptr->IntegrationElement(qr_quad_.Points()));
      // Fetch the transformation matrices for the gradients
      const Eigen::MatrixXd JinvT(
          geo_ptr->JacobianInverseGramian(qr_quad_.Points()));

      // Element matrix
      elem_mat_t mat(Nrsf_quad_, Nrsf_quad_);
      mat.setZero();

      // Loop over quadrature points
      for (int k = 0; k < Nqp_quad_; ++k) {
        // Evaluate diffusion and reaction coefficient
        // at quadrature points in actual cell
        const auto alphaval{alpha_(mapped_qpts.col(k))};
        const auto gammaval{gamma_(mapped_qpts.col(k))};

        const double w = qr_quad_.Weights()[k] * determinants[k];
        // Transformed gradients
        const auto trf_grad(JinvT.block(0, 2 * k, 2, 2) *
                            grad_quadpoint_quad_[k]);
        // Transformed gradients multiplied with coefficient
        const auto alpha_trf_grad(alphaval * trf_grad);
        mat += w * (alpha_trf_grad.transpose() * trf_grad +
                    (gammaval * rsf_quadpoints_quad_.col(k)) *
                        (rsf_quadpoints_quad_.col(k).transpose()));
      }
      return mat;
      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch
  return Eigen::MatrixXd(0, 0);
}

/**
 * @brief Local computation of general element vector for scalar finite elements
 *
 * @tparam SCALAR underlying scalar type, usually double or complex<double>
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<double(const Eigen::Vector2d &)>, which supplies
 *         the source function
 *
 * Computation is based on a quadrature rules supplied by the LehrFEM++
 * lf::quad::QuadRule module.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 */
template <typename SCALAR, typename FUNCTOR>
class ScalarFELocalLoadVector : public LocCompLagrFEPreprocessor {
 public:
  using elem_vec_t = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>;
  using ElemVec = const elem_vec_t;

  /** @brief Constructor, performs precomputations
   *
   * @param fe_trie local shape functions to be used on triangles
   * @param fe_quad local shape functions for quadrilaterals
   */
  ScalarFELocalLoadVector(const ScalarReferenceFiniteElement<SCALAR> &fe_tria,
                          const ScalarReferenceFiniteElement<SCALAR> &fe_quad,
                          FUNCTOR f);
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /*
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   *
   *
   */
  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  /** An object providing the source function */
  FUNCTOR f_;

 public:
  /*
   * @brief static variable for controlling (debugging) output
   */
  static unsigned int ctrl_;
  static const unsigned int kout_cell = 8;
  static const unsigned int kout_locvec = 16;
  static const unsigned int kout_dets = 32;
  static const unsigned int kout_loop = 64;
  static const unsigned int kout_qpts = 128;
};

template <typename SCALAR, typename FUNCTOR>
unsigned int ScalarFELocalLoadVector<SCALAR, FUNCTOR>::ctrl_ = 0;

// Constructors
template <typename SCALAR, typename FUNCTOR>
ScalarFELocalLoadVector<SCALAR, FUNCTOR>::ScalarFELocalLoadVector(
    const ScalarReferenceFiniteElement<SCALAR> &fe_tria,
    const ScalarReferenceFiniteElement<SCALAR> &fe_quad, FUNCTOR f)
    : LocCompLagrFEPreprocessor(fe_tria, fe_quad), f_(f) {}

template <typename SCALAR, typename FUNCTOR>
typename ScalarFELocalLoadVector<SCALAR, FUNCTOR>::ElemVec
ScalarFELocalLoadVector<SCALAR, FUNCTOR>::Eval(const lf::mesh::Entity &cell) {
  // Type for source function
  using source_fn_t = decltype(f_(Eigen::Vector2d::Zero()));
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
      // World coordinates of quadrature points
      const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr_tria_.Points()));
      SWITCHEDSTATEMENT(ctrl_, kout_qpts,
                        std::cout << "LOCVEC(Tria): Mapped quadrature points:\n"
                                  << mapped_qpts << std::endl);
      // Obtain the metric factors for the quadrature points
      const Eigen::VectorXd determinants(
          geo_ptr->IntegrationElement(qr_tria_.Points()));
      SWITCHEDSTATEMENT(ctrl_, kout_dets,
                        std::cout << "LOCVEC(Tria): Metric factors:\n"
                                  << determinants.transpose() << std::endl);
      // Element vector
      elem_vec_t vec(Nrsf_tria_);
      vec.setZero();

      // Loop over quadrature points
      for (int k = 0; k < Nqp_tria_; ++k) {
        // Source function value at quadrature point
        const auto fval = f_(mapped_qpts.col(k));
        SWITCHEDSTATEMENT(
            ctrl_, kout_loop,
            std::cout << "LOCVEC(Tria): ["
                      << qr_tria_.Points().col(k).transpose() << "] -> ["
                      << mapped_qpts.col(k).transpose() << "], f = " << fval
                      << ", weight = " << qr_tria_.Weights()[k] << std::endl);
        // Contribution of current quadrature point
        vec += (qr_tria_.Weights()[k] * determinants[k] * fval) *
               rsf_quadpoints_tria_.col(k);
      }
      SWITCHEDSTATEMENT(ctrl_, kout_locvec,
                        std::cout << "LOCVEC(Tria) = \n"
                                  << vec.transpose() << std::endl);
      return vec;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      // Quadrature points in world coordinates
      const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr_quad_.Points()));
      SWITCHEDSTATEMENT(ctrl_, kout_qpts,
                        std::cout << "LOCVEC(Quad): Mapped quadrature points:\n"
                                  << mapped_qpts << std::endl);
      // Obtain the metric factors for the quadrature points
      const Eigen::VectorXd determinants(
          geo_ptr->IntegrationElement(qr_quad_.Points()));
      SWITCHEDSTATEMENT(ctrl_, kout_dets,
                        std::cout << "LOCVEC(Quad): Metric factors:\n"
                                  << determinants.transpose() << std::endl);
      // Element vector
      elem_vec_t vec(Nrsf_quad_);
      vec.setZero();

      // Loop over quadrature points
      for (int k = 0; k < Nqp_quad_; ++k) {
        // Source function value at quadrature point
        const auto fval = f_(mapped_qpts.col(k));
        SWITCHEDSTATEMENT(
            ctrl_, kout_loop,
            std::cout << "LOCVEC(Quad): ["
                      << qr_quad_.Points().col(k).transpose() << "] -> ["
                      << mapped_qpts.col(k).transpose() << "], f = " << fval
                      << ", weight = " << qr_quad_.Weights()[k] << std::endl);
        // Contribution of current quadrature point
        vec += (qr_quad_.Weights()[k] * determinants[k] * fval) *
               rsf_quadpoints_quad_.col(k);
      }
      SWITCHEDSTATEMENT(ctrl_, kout_locvec,
                        std::cout << "LOCVEC(Quad) = \n"
                                  << vec.transpose() << std::endl);
      return vec;
      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch
  return Eigen::VectorXd(0);
}

}  // namespace lf::fe

#endif
