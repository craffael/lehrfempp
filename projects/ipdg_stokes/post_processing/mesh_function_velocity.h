#ifndef THESIS_POST_PROCESSING_MESH_FUNCTION_VELOCITY_H
#define THESIS_POST_PROCESSING_MESH_FUNCTION_VELOCITY_H

/**
 * @file mesh_function_velocity.h
 * @brief Compute the velocity from basis function coefficients of a vector
 * potential
 */

#include <lf/uscalfe/uscalfe.h>

namespace projects::ipdg_stokes {

namespace post_processing {

/**
 * @brief A MeshFunction returning the velocity computed from the basis function
 * coefficients of a vector potential
 * @tparam SCALAR_FE The scalar type used in the finite element space
 * @tparam SCALAR_COEFF The scalar type used in the coefficient vector
 */
template <typename SCALAR_FE, typename SCALAR_COEFF>
class MeshFunctionVelocity {
 public:
  /**
   * @brief Create a new MeshFunctionVelocity from a scalar finite element space
   * and a coefficient vector
   * @param fe_space A shared pointer to the scalar finite element space
   * containing the vector potential
   * @param dof_vector A vector containing the basis function coefficients of
   * the vector potential
   */
  MeshFunctionVelocity(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<SCALAR_FE>>
          fe_space,
      const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>& dof_vector)
      : grad_(fe_space, dof_vector) {}

  auto operator()(const lf::mesh::Entity& entity,
                  const Eigen::MatrixXd& local) const {
    // Compute the curl by rotating the gradient
    auto gradients = grad_(entity, local);
    Eigen::Matrix2d rot;
    rot << 0, 1, -1, 0;
    for (auto& g : gradients) {
      g = rot * g;
    }
    return gradients;
  }

 private:
  const lf::uscalfe::MeshFunctionGradFE<SCALAR_FE, SCALAR_COEFF> grad_;
};

template <typename T, typename SCALAR_COEFF>
MeshFunctionVelocity(std::shared_ptr<T>,
                     const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1>&)
    -> MeshFunctionVelocity<typename T::Scalar, SCALAR_COEFF>;

}  // end namespace post_processing

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_POST_PROCESSING_MESH_FUNCTION_VELOCITY_H
