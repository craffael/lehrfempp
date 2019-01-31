/**
 * @file
 * @brief simple functions demontrating use and capabilities of LehrFEM++ mesh module
 *        meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include "lecturedemomesh.h"

namespace lecturedemo {
using size_type = lf::base::size_type;
using dim_t = lf::base::dim_t;
using glb_idx_t = lf::base::glb_idx_t;
using sub_idx_t = lf::base::sub_idx_t;

/* SAM_LISTING_BEGIN_1 */
int traverseEntities(const lf::mesh::Mesh &mesh, dim_t codim) {
  LF_ASSERT_MSG((codim <= mesh.DimMesh()),
		"codim " << +codim << " too large");
  std::cout << "Mesh dimension = " << mesh.DimMesh()
            << ", iterating over entities of co-dim. " << codim << ": "
            << mesh.NumEntities(codim) << " exist" << std::endl;
  size_type cnt = 0;
  // Typical loop for running through all entities of a specific co-dimension
  for (const lf::mesh::Entity &entity : mesh.Entities(codim)) {
    // Print entity information including its \samemp{unique} index
    std::cout << cnt << ": Entity #" << mesh.Index(entity) << ": " << entity
              << std::endl;
    cnt++;
  }
  return cnt;
}
/* SAM_LISTING_END_1 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
std::pair<size_type, size_type> countCellTypes(const lf::mesh::Mesh &mesh) {
  size_type tria_cnt = 0, quad_cnt = 0;  // Counters
  // Loop over all cells (= co-dimension-0 entities) of the mesh
  for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
    // Fetch type information
    lf::base::RefEl ref_el{cell.RefEl()};
    // Test, if current cell is of a particular type
    switch (ref_el) {
      case lf::base::RefEl::kTria(): { tria_cnt++; break; }
      case lf::base::RefEl::kQuad(): { quad_cnt++; break; }
      default: { LF_VERIFY_MSG(false, "Unknown cell type"); }
    }}
  return {tria_cnt, quad_cnt};
}
/* SAM_LISTING_END_3 */
// clang-format on
// clang-format off
/* SAM_LISTING_BEGIN_2 */
void scanTopology(const lf::mesh::Mesh &mesh, dim_t codim) {
  LF_ASSERT_MSG((codim <= mesh.DimMesh()),
		"codim " << +codim << " too large");
  // loop over all entities of the specified codimension
  for (const lf::mesh::Entity &ent : mesh.Entities(codim)) {
    // Fetch topology type (TRIA or QUAD so far)
    const lf::base::RefEl ref_el{ent.RefEl()};
    // Print topological type and global index of the ent
    const glb_idx_t ent_idx = mesh.Index(ent);
    std::cout << ref_el << ": idx = " << ent_idx << std::endl;
    // Inspect sub-entities of any co-dimension
    for (dim_t sub_codim = 1; sub_codim <= mesh.DimMesh() - codim;
         ++sub_codim) {
      // Obtain iterator over sub-entities
      auto sub_ent_range = ent.SubEntities(sub_codim);
      size_type sub_cnt = 0;  // Counter for sub-entities
      // Loop over sub-entities, whose types and indices will be output
      for (const lf::mesh::Entity &subent : sub_ent_range) { // \Label[line]{st:1}
        std::cout << "\t rel. codim " << +sub_codim << " sub-ent "
		  << sub_cnt << ": " << subent << ", idx = "
		  << mesh.Index(subent) << std::endl;
        sub_cnt++;
      }}}}
/* SAM_LISTING_END_2 */
// clang-format on
// clang-format off
/* SAM_LISTING_BEGIN_4 */
void PrintGeometryInfo(const lf::mesh::Mesh &mesh, dim_t codim) {
  LF_ASSERT_MSG((codim <= mesh.DimMesh()),
		"codim " << +codim << " too large");
  // loop over all entities of the specified codimension
  for (const lf::mesh::Entity &ent : mesh.Entities(codim)) {
    // Number of nodes = number of corner points
    const size_type num_nodes = ent.RefEl().NumNodes();
    // Obtain pointer to geometry object associated with entity
    const lf::geometry::Geometry *geo_ptr = ent.Geometry();
    LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry!");
    // Fetch coordinates of corner points in packed format \cref{par:coords}
    Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
    LF_ASSERT_MSG(corners.rows() == geo_ptr->DimGlobal(),
                  "dimension mismatch for coordinate vectors");
    LF_ASSERT_MSG(corners.cols() == num_nodes, "#corners mismath");
    std::cout << ent.RefEl() << "(" << mesh.Index(ent) << ") pts: ";
    for (int l = 0; l < num_nodes; ++l) {
      std::cout << l << " =[" << corners.col(l).transpose() << "], ";
    }
    std::cout << std::endl;
  }}
/* SAM_LISTING_END_4 */
// clang-format on

/** @brief driver routine for LehrFEM++ demos for lecture
 */
void lecturedemomesh() {
  std::cout << "LehrFEM++ DEMO: mesh capabilities and functionality"
            << std::endl;
  // Set complete file path to the sample mesh in Gmsh format
  boost::filesystem::path here = __FILE__;
  auto mesh_file = here.parent_path() / "lecturedemomesh.msh";

  // clang-format off
  /* SAM_LISTING_BEGIN_5 */
  // Create a 2D mesh data structure from the information contained in the file
  // \texttt{mesh\_file}. A \com{factory object} is in charge of creating mesh
  // entities and has to be initialized first.
  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);// \Label[line]{lgm:1}
  lf::io::GmshReader reader(std::move(factory), mesh_file.string()); // \Label[line]{lgm:2}
  // Obtain pointer to read only mesh from the mesh reader object
  // Meshes in \lfpp are managed through \samemp{shared pointers}, see
  // \href{http://en.cppreference.com/w/cpp/memory/shared_ptr}{documentation}.
  std::shared_ptr<const lf::mesh::Mesh> mesh_ptr = reader.mesh();
  const lf::mesh::Mesh &mesh{*mesh_ptr};

  // Output general information on mesh; self-explanatory
  std::cout << "Mesh from file " << mesh_file.string() << ": ["
            << mesh.DimMesh() << ',' << mesh.DimWorld()
            << "] dim:" << std::endl;
  std::cout << mesh.NumEntities(0) << " cells, " << mesh.NumEntities(1)
            << " edges, " << mesh.NumEntities(2) << " nodes" << std::endl;
  /* SAM_LISTING_END_5 */
  // clang-format on

  // First demo: container functionalty of a mesh object
  for (lf::base::dim_t codim = 0; codim <= 2; ++codim) {
    (void)lecturedemo::traverseEntities(mesh, codim);
  }
  // Count number of instances of different cell types
  // NOLINTNEXTLINE
  auto [tria_cnt, quad_cnt] = lecturedemo::countCellTypes(mesh);
  std::cout << tria_cnt << " TRIA, " << quad_cnt << " QUAD cells" << std::endl;

  // Output information about topological relationships
  lecturedemo::scanTopology(mesh, 0);  // topolgy from the cell perspective
  lecturedemo::scanTopology(mesh, 1);  // topolgy from the edge perspective

  // Print information about the geometry of the mesh cells
  lecturedemo::PrintGeometryInfo(mesh, 0);
  // Print information about the geometry of the mesh edges
  lecturedemo::PrintGeometryInfo(mesh, 1);

  // Output rendering of mesh in LaTeX/TikZ format
  using lf::io::TikzOutputCtrl;
  lf::io::writeTikZ(
      mesh, "demomesh.tex",
      TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
          TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
          TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::ArrowTips);
}

}  // namespace lecturedemo
