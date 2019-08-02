/**
 * @file
 * @brief Implementation of a MeshFunction which represents a concrete
 * function from an approximation space.
 * @author Raffael Casagrande
 * @date   2019-01-12 06:05:36
 * @copyright MIT License
 */

#ifndef __4ee2d6e8004446558bc6d2186596e392
#define __4ee2d6e8004446558bc6d2186596e392

#include <memory>

#include "uniform_scalar_fe_space.h"

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @ingroup mesh_function
 * @brief A \ref mesh_function "MeshFunction" representing an element from a
 * ScalarUniformFESpace (e.g. solution of BVP)
 * @tparam SCALAR_FE The scalar type of the finite element basis functions.
 * @tparam SCALAR_COEFF The scalar type of the coefficient vector
 *
 * The MeshFunctionFE takes essentially two parameters:
 * - A ScalarUniformFESpace which defines the space of all approximation
 * functions.
 * - A Coefficient Vector which defines what element should be picked from the
 * FeSpace. This is a simple Eigen::Vector which defines the coefficient in
 * front of every basis function of the FeSpace (same ordering).
 *
 * @note A MeshFunctionFE can be evaluated on entities of all codimensions.
 *   The mesh function will simply use the ScalarReferenceFiniteElement from
 *   the feSpace.
 */
template <class SCALAR_FE, class SCALAR_COEFF>
class MeshFunctionFE {
 public:
  using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));

  /**
   * @brief Create a new MeshFunctionFE from a ScalarUniformFESpace and a
   * coefficient vector
   * @param fe_space the approximation space in which the function lies.
   * @param coeff_vector Defines the coefficients in front of the basis
   * functions of `fe_space`
   */
  MeshFunctionFE(
      std::shared_ptr<const UniformScalarFESpace<SCALAR_FE>> fe_space,
      const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& coeff_vector)
      : fe_space_(std::move(fe_space)), dof_vector_(coeff_vector) {
    for (auto& ref_el : {base::RefEl::kPoint(), base::RefEl::kSegment(),
                         base::RefEl::kTria(), base::RefEl::kQuad()}) {
      fe_[ref_el.Id()] = fe_space_->ShapeFunctionLayout(ref_el);
    }
  }

  /** Evaluate the mesh function on a MeshEntity */
  std::vector<Scalar> operator()(const lf::mesh::Entity& e,
                                 const Eigen::MatrixXd& local) const {
    auto sf_eval = fe_[e.RefEl().Id()]->EvalReferenceShapeFunctions(local);

    Eigen::Matrix<SCALAR_COEFF, 1, Eigen::Dynamic> local_dofs(1,
                                                              sf_eval.rows());
    auto global_dofs = fe_space_->LocGlobMap().GlobalDofIndices(e);
    for (Eigen::Index i = 0; i < sf_eval.rows(); ++i) {
      local_dofs(i) = dof_vector_(global_dofs[i]);
    }

    std::vector<Scalar> result(local.cols());
    Eigen::Map<Eigen::Matrix<Scalar, 1, Eigen::Dynamic>> temp(&result[0], 1,
                                                              local.cols());
    temp = local_dofs * sf_eval;
    return result;
  }

 private:
  std::shared_ptr<const UniformScalarFESpace<SCALAR_FE>> fe_space_;
  const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector_;
  std::array<std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR_FE>>, 5>
      fe_;
};

// deduction guide
template <class T, class SCALAR_COEFF>
MeshFunctionFE(std::shared_ptr<T>,
               const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>&)
    ->MeshFunctionFE<typename T::Scalar, SCALAR_COEFF>;

}  // namespace lf::uscalfe

#endif  // __4ee2d6e8004446558bc6d2186596e392
