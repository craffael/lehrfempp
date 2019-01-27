/**
 * @file
 * @brief simple functions demontrating use and capabilities of LehrFEM++
 * assemble module; meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include "lecturedemodof.h"

namespace lecturedemo {

// clang-format off
/* SAM_LISTING_BEGIN_1 */
void printDofInfo(const lf::assemble::DofHandler &dofh) {
  // Obtain pointer to the underlying mesh
  auto mesh = dofh.Mesh();
  // Number of degrees of freedom managed by the DofHandler object
  const lf::assemble::size_type N_dofs(dofh.NoDofs());
  std::cout << "DofHandler(" << dofh.NoDofs() << " dofs):" << std::endl;
  // Output information about dofs for entities of all co-dimensions
  for (lf::base::dim_t codim = 0; codim <= mesh->DimMesh(); codim++) {
    // Visit all entities of a codimension codim
    for (const lf::mesh::Entity &e : mesh->Entities(codim)) {
      // Fetch unique index of current entity supplied by mesh object
      const lf::base::glb_idx_t e_idx = mesh->Index(e);
      // Number of shape functions covering current entity
      const lf::assemble::size_type no_dofs(dofh.NoLocalDofs(e));
      // Obtain global indices of those shape functions ...
      lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> dofarray{
          dofh.GlobalDofIndices(e)};
      // and print them
      std::cout << e << ' ' << e_idx << ": " << no_dofs << " dofs = [";
      for (int loc_dof_idx = 0; loc_dof_idx < no_dofs; ++loc_dof_idx) {
        std::cout << dofarray[loc_dof_idx] << ' ';
      }
      std::cout << ']';
      // Also output indices of interior shape functions
      lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> intdofarray{
          dofh.InteriorGlobalDofIndices(e)};
      std::cout << " int = [";
      for (lf::assemble::gdof_idx_t int_dof : intdofarray) {
        std::cout << int_dof << ' ';
      }
      std::cout << ']' << std::endl;
    }
  }
  // List entities associated with the dofs managed by the current
  // DofHandler object
  for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
    const lf::mesh::Entity &e(dofh.Entity(dof_idx));
    std::cout << "dof " << dof_idx << " -> " << e << " " << mesh->Index(e)
              << std::endl;
  }
}  // end function printDofInfo

/* SAM_LISTING_END_1 */
// clang-format on

void lecturedemodof() {
  // Short name for 2d coordinate vectors
  using coord_t = Eigen::Vector2d;
  // Corner coordinates for quadrilateral
  using quad_coord_t = Eigen::Matrix<double, 2, 4>;
  // Corner coordinates for a triangle
  using tria_coord_t = Eigen::Matrix<double, 2, 3>;

  std::cout << "LehrFEM++ demo for dof handling and assembly facilities"
            << std::endl;
  // ======================================================================
  // I: Generate simple hybrid mesh comprising two cells, one quadrilateral,
  // one triangle

  // Obtain mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  mesh_factory_ptr->AddPoint(coord_t({0, 0}));      // point 0
  mesh_factory_ptr->AddPoint(coord_t({1, 0}));      // point 1
  mesh_factory_ptr->AddPoint(coord_t({0, 1}));      // point 2
  mesh_factory_ptr->AddPoint(coord_t({1, 1}));      // point 3
  mesh_factory_ptr->AddPoint(coord_t({1.5, 0.5}));  // point 4
  quad_coord_t quad_coord(2, 4);
  // First cell: the unit square
  quad_coord << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      lf::base::ForwardRange<const size_type>({0, 1, 3, 2}),
      std::make_unique<lf::geometry::Parallelogram>(quad_coord));
  // Second cell: an affine triangle
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({1, 3, 4}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Ready to build the mesh data structure
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  // Output mesh in LaTeX TikZ format
  using lf::io::TikzOutputCtrl;
  lf::io::writeTikZ(
      *mesh_p, "dofmesh.tex",
      TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
          TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
          TikzOutputCtrl::NodeNumbering);

  // ======================================================================
  // Step II: Initialize dof handler and divulge its information

  // Set number of shape functions associated with every entity type in order to
  // define a dof handler for a uniform layout of local shape functions.
  const size_type ndof_node = 1;  // One interior dof per node
  const size_type ndof_edge = 2;  // Two interior dof per edge
  const size_type ndof_tria = 1;  // One interior dof per triangle
  const size_type ndof_quad = 4;  // One interior dof per quad

  std::cout << "LehrFEM++ demo: assignment of global shape functions"
            << std::endl;
  std::cout << "#dof/vertex = " << ndof_node << std::endl;
  std::cout << "#dof/edge = " << ndof_edge << std::endl;
  std::cout << "#dof/triangle = " << ndof_tria << std::endl;
  std::cout << "#dof/quadrilateral = " << ndof_quad << std::endl;

  // Create a dof handler object describing a uniform distribution
  // of shape functions
  lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kPoint(), ndof_node},
               {lf::base::RefEl::kSegment(), ndof_edge},
               {lf::base::RefEl::kTria(), ndof_tria},
               {lf::base::RefEl::kQuad(), ndof_quad}});
  printDofInfo(dof_handler);
}

}  // namespace lecturedemo
