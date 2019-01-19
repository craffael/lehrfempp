/**
 * @file
 * @brief Implementation of a MeshFunction which represents a concrete
 * function from an approximation space.
 * @author Raffael Casagrande
 * @date   2019-01-12 06:05:36
 * @copyright MIT License
 */

#include <memory>

#include "scalar_uniform_fe_space.h"

namespace lf::uscalfe {

template <class SCALAR_FE, class SCALAR_COEFF>
class MeshFunctionFE {
 public:
  using Scalar = decltype(SCALAR_FE(0) * SCALAR_COEFF(0));

  MeshFunctionFE(
      const std::shared_ptr<ScalarUniformFESpace<SCALAR_FE>>& fe_space,
      const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector)
      : fe_space_(fe_space), dof_vector_(dof_vector) {}

  std::vector<Scalar> operator()(const lf::mesh::Entity& e,
                                 const Eigen::MatrixXd& local) const {
    auto shape_functions = fe_space_->ShapeFunctionLayout(e.RefEl());
    auto sf_eval = shape_functions->EvalReferenceShapeFunctions(local);

    Eigen::Matrix<SCALAR_COEFF, 1, Eigen::Dynamic> local_dofs(1,
                                                              sf_eval.rows());
    auto global_dofs = fe_space_->LocGlobMap().GlobalDofIndices(e);
    for (int i = 0; i < sf_eval.rows(); ++i) {
      local_dofs(i) = dof_vector_(global_dofs[i]);
    }

    std::vector<Scalar> result(local.cols());
    Eigen::Map<Eigen::Matrix<Scalar, 1, Eigen::Dynamic>> temp(&result[0], 1,
                                                              local.cols());
    temp = local_dofs * sf_eval;
    return result;
  }

 private:
  std::shared_ptr<ScalarUniformFESpace<SCALAR_FE>> fe_space_;
  const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector_;
};

// deduction guide
template <class T, class SCALAR_COEFF>
MeshFunctionFE(std::shared_ptr<T>,
               const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>&)
    ->MeshFunctionFE<typename T::Scalar, SCALAR_COEFF>;

}  // namespace lf::uscalfe
