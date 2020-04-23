#ifndef PROJECTS_DPG_SUB_ELEMENT_VECTOR_PROVIDER
#define PROJECTS_DPG_SUB_ELEMENT_VECTOR_PROVIDER

/**
 * @file
 * @brief Interface Class providing sub element vectors associated
 * with linear forms of a component of a cartesian/product space.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"

namespace projects::dpg {

/**
 * @brief Interface class providing element vectors associated
 * with linear forms of a component of a cartesian/product space.
 *
 * @tparam SCALAR type for the entries of the elemnt matrices. Ususlly 'double'
 * @note This interface class complies with the requirements for the template
 * parameter `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 *
 * This class provides sub element vectors used to construct element vectors
 * associated to linear forms \f$ l: V \rightarrow \mathbb R \f$ on a product
 * space (c.f. ProductUniformFEDofHandler, ProductUniformFESpace)
 *
 *
 * \f[ V = V_0 \times V_1 \times \dots \times V_{m-1} \f]
 *
 * that have the following structure
 *
 * \f[ l((v_1, \dots v_{m-1})) = \sum_{k} l_k(v_{j_k}) \f]
 *
 * with
 * \f[ l_k :  V_{j_k} \rightarrow \mathbb{R} \f]
 *
 * This class allows the evaluation of element vectors associated to such
 * "building block" linear forms \f$ l_k \f$ and provides information about the
 * occuring component \f$ v_{j_k} \f$. Using this information, the evaluated
 * element vectors are used as sub vectors to construct the element matrix
 * associated to the linear form \f$ l \f$ in the class
 * ProductElementVectorProvider.
 *
 */
template <typename SCALAR>
class SubElementVectorProvider {
 public:
  using size_type = lf::uscalfe::size_type;
  /** @brief internal type for element vectors*/
  using elem_vec_t = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>;
  /** Return type of the @ref Eval() method */
  using ElemVec = elem_vec_t;

  SubElementVectorProvider() = default;
  virtual ~SubElementVectorProvider() = default;

  SubElementVectorProvider(const SubElementVectorProvider&) = delete;
  SubElementVectorProvider(SubElementVectorProvider&&) noexcept = default;
  SubElementVectorProvider& operator=(const SubElementVectorProvider&) = delete;
  SubElementVectorProvider& operator=(SubElementVectorProvider&&) = delete;

  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }

  /**
   * @brief main routine for the computation of (sub) element vectors
   * @param cell refernce to the cell for which the (sub) element
   * vector should be computed
   * @return a small vector, containing the "building block" element vector.
   */
  virtual ElemVec Eval(const lf::mesh::Entity& cell) = 0;

  /**
   * @brief returns the index of the test space component \f$ v_{j_k} \f$ which
   * is the test space for the linear form \f$ l_k \f$
   */
  [[nodiscard]] virtual size_type TestComponent() const = 0;
};

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_SUB_ELEMENT_VECTOR_PROVIDER
