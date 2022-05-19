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
 * @brief Computes the squared L2Norm of some MeshFunction (type MF) f
 *
 * @tparam MF type of the mesh function to take to norm form
 * @tparam SQ_F type of the square function sq_f
 * @param mesh_p pointer to the mesh
 * @param f meshfunction to evaluate
 * @param sq_f square function for the output of f (sq_f(x) -> |x|^2)
 * @param quadrule quadrature rule used for each cell in the mesh
 *
 * @returns pair with the squared norm of the inputfunction f on the mesh and
 *  CodimMeshDataSet containing the cellwise integrals over the squared errors
 *
 * Iterates over all cells in the mesh and compute the
 * squared integral locally, then add up and square
 *
 */
template <typename MF, typename SQ_F>
std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2norm(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p, const MF& f,
    const SQ_F& sq_f, const lf::quad::QuadRule& quadrule) {
  // store the intermediate squared sums of cells
  double squared_sum = 0;

  // create codim mesh dataset to store the cellwise integrals
  lf::mesh::utils::CodimMeshDataSet<double> cellErrors(mesh_p, 0);

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
    double cell_error = 0;
    for (lf::base::size_type i = 0; i < n_points; i++) {
      cell_error += det[i] * weights[i] * sq_f(values[i]);
    }
    cellErrors(*e) = cell_error;
    squared_sum += cell_error;
  }

  return std::make_pair(squared_sum, cellErrors);
}
}  // namespace projects::hldo_sphere::post_processing
#endif  // THESIS_NORMS_L2NORM_H
