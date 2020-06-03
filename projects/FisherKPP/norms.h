/** @file norms.h
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 05.05.20
 *  @copyright ETH Zurich
 */

#include <Eigen/Core>

#include <lf/mesh/mesh.h>

namespace FisherKPP {

double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p);

Eigen::VectorXd reduce(const Eigen::VectorXd &mu, unsigned int N);

} /* namespace FisherKPP */
