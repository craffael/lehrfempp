#ifndef PROJECTS_DPG_SUB_ELEMENT_MATRIX_PROVIDER
#define PROJECTS_DPG_SUB_ELEMENT_MATRIX_PROVIDER

/**
 * @file
 * @brief Interface Class providing sub element matrices associated
 * with bilinear forms between components of cartesian/product spaces.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"

namespace projects::dpg {

/**
 * @brief Interface class providing sub element matrices associated
 * with bilinear forms between components of Cartesian product spaces.
 *
 * @tparam SCALAR type for the entries of the element matrices. Ususally
 * 'double'
 * @note This interface class complies with the type requirements for the
 * template argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * This class provides sub element matrices used to construct element matrices
 * associated to bilinear forms \f$ b: U \times V \rightarrow \mathbb{R} \f$
 * between two product spaces (c.f. ProductUniformFEDofHandler,
 * ProductUniformFESpace)
 *
 * \f[ U = U_0 \times U_1 \times \dots \times U_{n-1} \f]
 * \f[ V = V_0 \times V_1 \times \dots \times V_{m-1} \f]
 *
 * that have the following structure
 *
 * \f[ b((u_1, \dots, u_{n-1}),(v_1, \dots v_{m-1})) = \sum_{k}
 * b_k(u_{i_k},v_{j_k}) \f]
 *
 * with
 * \f[ b_k : U_{i_k} \times V_{j_k} \rightarrow \mathbb{R} \f]
 *
 * This class allows the evaluation of element matries associated to such
 * "building block" bilinear forms \f$ b_k \f$ and provides information about
 * the occuring components \f$ u_{i_k}, v_{j_k} \f$. Using this information, the
 * evaluated element matrices are used as sub matrices to construct the element
 * matrices associated to the bilinear form \f$ b \f$ in the class
 * ProductElementMatrixProvider.
 *
 */
template <typename SCALAR>
class SubElementMatrixProvider {
 public:
  using size_type = lf::uscalfe::size_type;
  /** @brief internal type for element matrices*/
  using elem_mat_t = Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>;
  /** Return type of the @ref Eval() method */
  using ElemMat = elem_mat_t;

  SubElementMatrixProvider() = default;
  virtual ~SubElementMatrixProvider() = default;

  /** @brief standard constructors */
  SubElementMatrixProvider(const SubElementMatrixProvider&) = delete;
  SubElementMatrixProvider(SubElementMatrixProvider&&) noexcept = default;
  SubElementMatrixProvider& operator=(const SubElementMatrixProvider&) = delete;
  SubElementMatrixProvider& operator=(SubElementMatrixProvider&&) = delete;

  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief main routine for the computation of (sub) element
   * matrices
   * @param cell refernce to the cell for which the  element
   * matrix should be computed
   * @return a small dense, containing the element matrix.
   */
  virtual ElemMat Eval(const lf::mesh::Entity& cell) = 0;

  /**
   * @brief returns the index of the trial space component \f$ u_{i_k} \f$ which
   * is the trial space for the bilinear form \f$ b_k \f$
   */
  [[nodiscard]] virtual size_type TrialComponent() const = 0;

  /**
   * @brief returns the index of the test space component \f$ v_{j_k} \f$ which
   * is the test space for the bilinear form \f$ b_k \f$
   */
  [[nodiscard]] virtual size_type TestComponent() const = 0;
};

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_SUB_ELEMENT_MATRIX_PROVIDER
