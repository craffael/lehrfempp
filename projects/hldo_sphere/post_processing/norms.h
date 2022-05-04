#ifndef THESIS_NORMS_L2NORM_H
#define THESIS_NORMS_L2NORM_H

/**
 * @file norms.h
 *
 * @brief provides norms for the evaluation of the errors
 */

namespace projects::hldo_sphere::post_processing {
/**
 *
 * @brief Computes the L2Norm of some MeshFunction (type MF) f
 *
 * @tparam MF type of the mesh function to take to norm form
 * @tparam SQ_F type of the square function sq_f
 * @param mesh_p pointer to the mesh
 * @param f meshfunction to evaluate
 * @param sq_f square function for the output of f (sq_f(x) -> |x|^2)
 * @param quadrule quadrature rule used for each cell in the mesh
 *
 * @returns the squared norm of the inputfunction f on the mesh
 *
 * Iterates over all cells in the mesh and compute the
 * squared integral locally, then add up and square
 *
 */
template <typename MF, typename SQ_F>
double L2norm(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p, const MF& f,
              const SQ_F& sq_f, const lf::quad::QuadRule& quadrule) {
  // store the intermediate squared sums of cells
  double squared_sum = 0;

  // get ref quad points and weights for each cell
  const Eigen::MatrixXd loc_points = quadrule.Points();
  const Eigen::VectorXd weights = quadrule.Weights();
  const lf::base::size_type n_points = quadrule.NumPoints();

  // loop over all cells in the mesh
  for (const lf::mesh::Entity* e : mesh_p->Entities(0)) {
    // get determinante of the pullpack this is constant for all points
    Eigen::VectorXd det{e->Geometry()->IntegrationElement(loc_points)};
    // get funciton values
    auto values = f(*e, loc_points);
    // compute local quadrature of the squared norm
    for (lf::base::size_type i = 0; i < n_points; i++) {
      squared_sum += det[i] * weights[i] * sq_f(values[i]);
    }
  }

  return std::sqrt(squared_sum);
}
}  // namespace projects::hldo_sphere::post_processing
#endif  // THESIS_NORMS_L2NORM_H
