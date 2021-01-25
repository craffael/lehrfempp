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

#include "scalar_fe_space.h"

namespace lf::fe {

/**
 * @headerfile lf/fe/fe.h
 * @ingroup mesh_function
 * @brief A \ref mesh_function "MeshFunction" representing an element from a
 * ScalarUniformFESpace (e.g. solution of BVP)
 * @tparam SCALAR_FE The scalar type of the finite element basis functions.
 * @tparam SCALAR_COEFF The scalar type of the coefficient vector
 *
 * The MeshFunctionFE takes essentially two parameters:
 * - A ScalarFESpace which defines the space of all approximation
 * functions.
 * - A Coefficient Vector which defines what element should be picked from the
 * FeSpace. This is a simple Eigen::Vector which defines the coefficient in
 * front of every basis function of the FeSpace (same ordering).
 *
 * @note A MeshFunctionFE can be evaluated on entities of all codimensions.
 *   The mesh function will simply use the ScalarReferenceFiniteElement from
 *   the feSpace.
 *
 * ### Use case
 *
 * @snippet mesh_function_binary.cc mffedemo
 *
 */
template <class SCALAR_FE, class SCALAR_COEFF>
class MeshFunctionFE {
 public:
  // Why this? Because we can use real-valued finite element spaces to represent
  // complex-valued functions by using complex degrees of freedom!
  using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));

  /**
   * @brief Create a new MeshFunctionFE from a @ref ScalarFESpace and a
   * coefficient vector
   * @param fe_space the approximation space in which the function lies.
   * @param coeff_vector Defines the coefficients in front of the basis
   * functions of `fe_space`
   *
   * @warning This constructor just takes a reference to the vector of basis
   * expansion coefficients. Thus, this vector has to be "kept alive" as long as
   * the mesh function exists.
   */
  MeshFunctionFE(
      std::shared_ptr<const ScalarFESpace<SCALAR_FE>> fe_space,
      const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& coeff_vector)
      : fe_space_(std::move(fe_space)), dof_vector_(coeff_vector) {}

  /** @brief Evaluate the mesh function on a MeshEntity
   *
   * @param e the relevant mesh entity
   * @param local _reference coordinates_ of evalation points passed in the
   * columns of a matrix
   * @return values of the function at evaluation points
   */
  std::vector<Scalar> operator()(const lf::mesh::Entity& e,
                                 const Eigen::MatrixXd& local) const {
    // Obtain values of all shape functions in the evaluation points
    auto sf_eval =
        fe_space_->ShapeFunctionLayout(e)->EvalReferenceShapeFunctions(local);
    // Extract the d.o.f.s for the current entity from the vector of global
    // d.o.f. values
    Eigen::Matrix<SCALAR_COEFF, 1, Eigen::Dynamic> local_dofs(1,
                                                              sf_eval.rows());
    auto global_dofs = fe_space_->LocGlobMap().GlobalDofIndices(e);
    for (Eigen::Index i = 0; i < sf_eval.rows(); ++i) {
      local_dofs(i) = dof_vector_(global_dofs[i]);
    }
    // Trick to combine Eigen data types with STL containers
    std::vector<Scalar> result(local.cols());
    Eigen::Map<Eigen::Matrix<Scalar, 1, Eigen::Dynamic>> temp(&result[0], 1,
                                                              local.cols());
    temp = local_dofs * sf_eval;
    return result;
  }

  /**
   * @brief Convenience method to retrieve the underlying mesh
   * @returns The mesh on which this mesh function is defined.
   */
  [[nodiscard]] std::shared_ptr<const mesh::Mesh> getMesh() const {
    return fe_space_->Mesh();
  }

  /**
   * @brief Convenience method to retrieve the finite element space in which the
   * mesh function lives
   * @returns shared_ptr to ScalarFESpace in which the mesh function
   * lives.
   */
  [[nodiscard]] std::shared_ptr<const ScalarFESpace<SCALAR_FE>> getFESpace()
      const {
    return fe_space_;
  }

 private:
  std::shared_ptr<const ScalarFESpace<SCALAR_FE>> fe_space_;
  const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector_;
};

// deduction guide
template <class T, class SCALAR_COEFF>
MeshFunctionFE(std::shared_ptr<T>,
               const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>&)
    ->MeshFunctionFE<typename T::Scalar, SCALAR_COEFF>;

}  // namespace lf::fe

#endif  // __4ee2d6e8004446558bc6d2186596e392
