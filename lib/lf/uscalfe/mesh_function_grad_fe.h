/**
 * @file
 * @brief Represents the gradient of an element of a finite element space.
 * @author Raffael Casagrande
 * @date   2019-01-19 06:39:07
 * @copyright MIT License
 */

#ifndef __b6997524e2834b5b8e4bba019fb35cc6
#define __b6997524e2834b5b8e4bba019fb35cc6
#include "scalar_uniform_fe_space.h"

namespace lf::uscalfe {

template <class SCALAR_FE, class SCALAR_COEFF>
class MeshFunctionGradFE {
 public:
  using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));

  MeshFunctionGradFE(
      const std::shared_ptr<ScalarUniformFESpace<SCALAR_FE>>& fe_space,
      const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector)
      : fe_space_(fe_space), dof_vector_(dof_vector) {}

  std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> operator()(
      const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const {
    auto shape_functions = fe_space_->ShapeFunctionLayout(e.RefEl());
    auto grad_sf_eval =
        shape_functions->GradientsReferenceShapeFunctions(local);

    Eigen::Matrix<SCALAR_COEFF, 1, Eigen::Dynamic> local_dofs(
        1, grad_sf_eval.rows());
    auto global_dofs = fe_space_->LocGlobMap().GlobalDofIndices(e);
    for (int i = 0; i < grad_sf_eval.rows(); ++i) {
      local_dofs(i) = dof_vector_(global_dofs[i]);
    }

    // gradients w.r.t. reference element coordinates
    auto local_grads = (local_dofs * grad_sf_eval).eval();

    auto jac_t = e.Geometry()->JacobianInverseGramian(local);
    auto dim_local = e.RefEl().Dimension();
    std::vector<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> result(local.cols());
    for (int i = 0; i < result.size(); ++i) {
      result[i] = jac_t.block(0, dim_local * i, jac_t.rows(), dim_local) *
                  local_grads.block(0, i * dim_local, 1, dim_local).transpose();
    }

    return result;
  }

 private:
  std::shared_ptr<ScalarUniformFESpace<SCALAR_FE>> fe_space_;
  const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector_;
};

// deduction guide
template <class T, class SCALAR_COEFF>
MeshFunctionGradFE(std::shared_ptr<T>,
                   const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>&)
    ->MeshFunctionGradFE<typename T::Scalar, SCALAR_COEFF>;

}  // namespace lf::uscalfe

#endif  // __b6997524e2834b5b8e4bba019fb35cc6
