#ifndef THESIS_MESH_FUNCTION_WHITNEY_TWO_H
#define THESIS_MESH_FUNCTION_WHITNEY_TWO_H

#include <lf/uscalfe/uscalfe.h>

/**
 * @file mesh_function_whitney_two.h
 *
 * @brief Provides utilities for the evaluation of errors
 */

namespace projects::hldo_sphere::post_processing {

class MeshFunctionWhitneyTwo {
 public:
  /**
   * @brief basic constructor
   *
   * @param mu vector containing the basis function expansion coefficiants in
   * global ordering
   * @param mesh containing the mesh on which to evaluate
   */
  MeshFunctionWhitneyTwo(const Eigen::VectorXd& mu,
                         const std::shared_ptr<const lf::mesh::Mesh> mesh)
      : mu_(mu), mesh_(mesh) {
    // make sure mu has the right size
    lf::base::size_type n = mesh->NumEntities(0);
    int mu_size = mu.size();
    LF_VERIFY_MSG(
        n == mu_size,
        "Not the right number of basis expansion coefficiants in mu expected: "
            << n << " given: " << mu_size);
  };

  /**
   *
   * @brief Evaluates the whitney one form basis function at reference points on
   * a given cell
   *
   * @param e on which to evaluate
   * @param local coordinates on the reference triangle for the
   * evaluation. Each column contains one evaluation point
   *
   *
   * @returns a vector of scalars containing the mu value
   * for the given mesh since the function is piecewise constant.
   *
   * And returns a vector containing the values as follows
   *
   * @f[
   *   \mu_{e}
   * ]@f
   *
   * where the $\lambda$ are the barycentric basis functions
   *
   *
   * @note only triangles are supported
   *
   */
  std::vector<double> operator()(const lf::mesh::Entity& e,
                                 const Eigen::MatrixXd& local) {
    // Only triangles are supported
    LF_VERIFY_MSG(e.RefEl() == lf::base::RefEl::kTria(),
                  "Unsupported cell type " << e.RefEl());

    lf::base::size_type global_index = mesh_->Index(e);

    std::vector<double> vals(local.cols());
    for (int i = 0; i < local.cols(); i++) {
      vals[i] = mu_(global_index);
    }

    return vals;
  }

 private:
  const Eigen::VectorXd& mu_;
  const std::shared_ptr<const lf::mesh::Mesh> mesh_;
};

}  // namespace projects::hldo_sphere::post_processing
#endif  // THESIS_MESH_FUNCTION_WHITNEY_TWO_H
