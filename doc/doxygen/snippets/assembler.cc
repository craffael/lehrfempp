#include <lf/assemble/assemble.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf {

void matrix() {
  //! [matrix_usage]
  // initialize triangular 2d mesh somehow
  std::shared_ptr<mesh::Mesh> mesh;

  // create a dof handler which assigns one global dof to every node of the mesh
  assemble::UniformFEDofHandler dofh(mesh, {{base::RefEl::kPoint(), 1}});

  // initialize EntityMatrixProvider that defines the local Laplace element
  // matrix (only for triangular meshes):
  uscalfe::LinearFELaplaceElementMatrix entity_matrix_provider;

  // setup a COOMatrix into which the global matrix will be assembled:
  assemble::COOMatrix<double> lhs(dofh.NoDofs(), dofh.NoDofs());

  // assemble the global Laplace matrix (iterate over all entities with
  // codim=2):
  assemble::AssembleMatrixLocally<assemble::COOMatrix<double>>(
      2, dofh, dofh, entity_matrix_provider, lhs);
  //![matrix_usage]
}

void vector() {
  //![vector_usage]
  // initialize a 2d mesh somehow
  std::shared_ptr<mesh::Mesh> mesh;

  // define a EntityVectorProvider which returns the volume for every mesh cell
  struct VolumeEntityVectorProvider {
    bool isActive(const lf::mesh::Entity& e) const { return e.Codim() == 0; }

    Eigen::VectorXd Eval(const lf::mesh::Entity& e) const {
      return Eigen::VectorXd::Constant(1, geometry::Volume(*e.Geometry()));
    }
  } entity_vector_provider;

  // define a dof handler which assigns a dof to every mesh cell
  assemble::UniformFEDofHandler dofh(
      mesh, {{base::RefEl::kTria(), 1}, {base::RefEl::kQuad(), 1}});

  // initialize the output vector:
  Eigen::VectorXd x(dofh.NoDofs());

  // assemble the global vector over entities with codim=0:
  AssembleVectorLocally(0, dofh, entity_vector_provider, x);

  // now x[dofh.GlobalDofIndices(e)[0]] will contain the volume of mesh cell e
  //![vector_usage]
}

}  // namespace lf
