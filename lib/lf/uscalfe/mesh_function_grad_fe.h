/**
 * @file
 * @brief Represents the gradient of an element of a finite element space.
 * @author Raffael Casagrande
 * @date   2019-01-19 06:39:07
 * @copyright MIT License
 */

#ifndef __b6997524e2834b5b8e4bba019fb35cc6
#define __b6997524e2834b5b8e4bba019fb35cc6
#include "uniform_scalar_fe_space.h"

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @ingroup mesh_function
 * @brief A \ref mesh_function "MeshFunction" representing the gradient of a
 * function from a finite element space (e.g. gradient of a solution of BVP).
 * @tparam SCALAR_FE The scalar type of the finite element basis functions.
 * @tparam SCALAR_COEFF The scalar type of the coefficient vector
 *
 * The MeshFunctionGradFE takes essentially two parameters:
 * - A ScalarUniformFESpace which defines the space of all approximation
 * functions. This should be a Lagrnagian finite element space of globally
 * continuous functions. Otherwise the evaluation for non-cell entities is not
 * possible.
 * - A Coefficient Vector which defines what element should be picked from the
 * FeSpace. This is a simple Eigen::Vector which defines the coefficient in
 * front of every basis function of the FeSpace (same ordering).
 *
 * @note A MeshFunctionGradFE can be evaluated on entities of all codimensions
 * except for points.
 *
 * @note The gradient that is returned by this mesh function is w.r.t. to the
 * global coordinate system of the mesh, i.e. the gradients of the shape
 * function are multiplied by Geometry::JacobianInverseGramian(). Hence, what is
 * returned is a "tangential gradient". For instance, if the evaluation operator
 * is invoked for an edge of a 2D mesh, the gradient of the finite element
 * function in the direction of that edge is returned. If the evaluation is
 * requested for a (curved) triangle in 3D, then a gradient with 3 components is
 * retrurned, which represents the tangential gradient of the finite element
 * function _in global 3D Euclidean coordinates_.
 */
template <class SCALAR_FE, class SCALAR_COEFF>
class MeshFunctionGradFE {
 public:
  using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));

  /**
   * @brief Create a new MeshFunctionGradFE from a ScalarUniformFESpace and a
   * coefficient vector
   * @param fe_space the approximation space in which the function lies.
   * @param dof_vector Defines the coefficients in front of the basis
   * functions of `fe_space`
   */
  MeshFunctionGradFE(
      std::shared_ptr<const UniformScalarFESpace<SCALAR_FE>> fe_space,
      const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector)
      : fe_space_(std::move(fe_space)), dof_vector_(dof_vector) {
    // Loop over entity types for which evaluations are supported. Then collect
    // information about local shape functions associated with those entity
    // types.
    for (auto& ref_el : {base::RefEl::kSegment(), base::RefEl::kTria(),
                         base::RefEl::kQuad()}) {
      fe_[ref_el.Id()] = fe_space_->ShapeFunctionLayout(ref_el);
    }
  }

  /**
   * @brief Convenience method to retrieve the underlying mesh
   * @return The mesh on which this mesh function is defined.
   */
  [[nodiscard]] std::shared_ptr<const mesh::Mesh> getMesh() const {
    return fe_space_->Mesh();
  }

  /**
   * @brief Convenience method to retrieve the finite element space to which the
   * original function belongs (i.e. before taking the gradient)
   * @return shared_ptr to UniformScalarFESpace to which the original function
   * belongs (i.e. before taking the gradient)
   */
  [[nodiscard]] std::shared_ptr<const UniformScalarFESpace<SCALAR_FE>>
  getFESpace() const {
    return fe_space_;
  }

  /** @brief Evaluate the mesh function on a MeshEntity
   *
   * @param e reference to a mesh entity
   * @param local matrix whose columns contain reference coordinates of the
   *          evaluation points on the entity e
   * @return array of column vectors containing the gradients.
   *
   */
  std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> operator()(
      const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const {
    // Access information on local shape functions for the entity
    const unsigned int type_id = e.RefEl().Id();
    LF_ASSERT_MSG(fe_[type_id] != nullptr, "Missing LSF information for " << e);
    auto grad_sf_eval = fe_[type_id]->GradientsReferenceShapeFunctions(local);
    // Fetch the coefficients of the global shape functions associated with the
    // current mesh entity
    Eigen::Matrix<SCALAR_COEFF, 1, Eigen::Dynamic> local_dofs(
        1, grad_sf_eval.rows());
    auto global_dofs = fe_space_->LocGlobMap().GlobalDofIndices(e);
    for (Eigen::Index i = 0; i < grad_sf_eval.rows(); ++i) {
      local_dofs(i) = dof_vector_(global_dofs[i]);
    }
    // gradients w.r.t. reference element coordinates \hat{x}
    auto local_grads = (local_dofs * grad_sf_eval).eval();
    // Transform to Cartesian coordinates 
    auto jac_t = e.Geometry()->JacobianInverseGramian(local);
    auto dim_local = e.RefEl().Dimension();
    std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> result(local.cols());
    // Transform all the local gradients in the evaluation points 
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = jac_t.block(0, dim_local * i, jac_t.rows(), dim_local) *
                  local_grads.block(0, i * dim_local, 1, dim_local).transpose();
    }
    return result;
  }

 private:
  /** @brief Pointer to underlying Lagrangian finite element space */
  std::shared_ptr<const UniformScalarFESpace<SCALAR_FE>> fe_space_;
  /** @brief Basis expansion coefficient vector for finite-element function */
  const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector_;
  /** @brief Pointers to descriptions of local shape functions */
  std::array<std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR_FE>>, 5>
      fe_;
};

// deduction guide
template <class T, class SCALAR_COEFF>
MeshFunctionGradFE(std::shared_ptr<T>,
                   const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>&)
    ->MeshFunctionGradFE<typename T::Scalar, SCALAR_COEFF>;

}  // namespace lf::uscalfe

#endif  // __b6997524e2834b5b8e4bba019fb35cc6
