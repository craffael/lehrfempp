#ifndef THESIS_MESH_FUNCTION_WHITNEY_ZERO_H
#define THESIS_MESH_FUNCTION_WHITNEY_ZERO_H

#include <lf/uscalfe/uscalfe.h>

/**
 * @file mesh_function_whitney_zero.h
 *
 * @brief Provides utilities for the evaluation of errors
 */

namespace projects::hldo_sphere::post_processing {

class MeshFunctionWhitneyZero {
 public:
  /**
   * @brief basic constructor
   *
   * @param mu vector containing the basis function expansion coefficiants in
   * global ordering
   * @param mesh containing the mesh on which to evaluate
   */
  MeshFunctionWhitneyZero(const Eigen::VectorXd& mu,
                          const std::shared_ptr<const lf::mesh::Mesh> mesh)
      : mu_(mu), mesh_(mesh) {
    // make sure mu has the right size
    lf::base::size_type n = mesh->NumEntities(2);
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
   * @returns a vector of doubles containing the results for the
   * evaluation of the linear combination of barycentric basis functions at the
   * input coordintes in `local`.
   *
   *
   * And returns a vector containing the values as follows
   *
   * @f[
   *          \mu_{l_0} \lambda_{0}(\bm{x}_i) + \mu_{l_1}\lambda_1(\bm{x}_i) +
   *              \mu_{l_2}\lambda_2(\bm{x}_i)
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

    // get the lambdas
    const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;

    // get lambda values
    const Eigen::MatrixXd lambda = hat_func.EvalReferenceShapeFunctions(local);

    // define return vector
    std::vector<double> vals(local.cols());

    // get global indices of the edges
    std::vector<lf::base::size_type> global_indices(3);
    auto points = e.SubEntities(2);
    for (lf::base::size_type i = 0; i < 3; i++) {
      global_indices[i] = mesh_->Index(*points[i]);
    }

    // for each point evaluate the basis functions and add the
    // values together
    for (int i = 0; i < local.cols(); i++) {
      vals[i] = mu_(global_indices[0]) * lambda(0) +
                mu_(global_indices[1]) * lambda(1) +
                mu_(global_indices[2]) * lambda(2);
    }

    return vals;
  }

 private:
  const Eigen::VectorXd& mu_;
  const std::shared_ptr<const lf::mesh::Mesh> mesh_;
};

}  // namespace projects::hldo_sphere::post_processing
#endif  // THESIS_MESH_FUNCTION_WHITNEY_ONE_H
