/**
 * @file
 * @brief Implementation of PrintInfo functions
 * @author Raffael Casagrande
 * @date   2018-07-01 01:33:21
 * @copyright MIT License
 */

#include "print_info.h"

namespace lf::mesh::utils {
void PrintInfo(const Mesh &mesh, std::ostream &o) {
  using dim_t = Mesh::dim_t;
  using size_type = Mesh::size_type;

  const dim_t dim_mesh = mesh.DimMesh();
  const dim_t dim_world = mesh.DimWorld();
  o << "Mesh of dimension " << dim_mesh << ", ambient dimension " << dim_world
    << std::endl;

  // Lopp over codimensions
  for (int co_dim = dim_mesh; co_dim >= 0; co_dim--) {
    const size_type no_ent = mesh.Size(co_dim);
    o << "Co-dimension " << co_dim << ": " << no_ent << " entities"
      << std::endl;

    // Loop over entities
    for (const Entity &e : mesh.Entities(co_dim)) {
      size_type e_idx = mesh.Index(e);
      dim_t e_codim = e.Codim();
      const geometry::Geometry *e_geo_ptr = e.Geometry();
      lf::base::RefEl e_refel = e.RefEl();

      LF_VERIFY_MSG(e_geo_ptr,
                    co_dim << "-entity " << e_idx << ": missing geometry");
      LF_VERIFY_MSG(e_geo_ptr->DimLocal() == dim_mesh - co_dim,
                    co_dim << "-entity " << e_idx << ": wrong dimension");
      LF_VERIFY_MSG(e_geo_ptr->RefEl() == e_refel,
                    co_dim << "-entity " << e_idx << ": refEl mismatch");
      LF_VERIFY_MSG(e_codim == co_dim,
                    co_dim << "-entity " << e_idx << " co-dimension mismatch");
      const Eigen::MatrixXd &ref_el_corners(e_refel.NodeCoords());
      o << "entity " << e_idx << " (" << e_refel << "): ";

      // Loop over local co-dimensions
      for (int l = 1; l <= dim_mesh - co_dim; l++) {
        o << "rel codim-" << l << " subent: [";
        // Fetch subentities of co-dimension l
        base::RandomAccessRange<const Entity> sub_ent_range(e.SubEntities(l));
        for (const Entity &sub_ent : sub_ent_range) {
          o << mesh.Index(sub_ent) << ' ';
        }
        o << "], ";
      }
      o << std::endl << e_geo_ptr->Global(ref_el_corners) << std::endl;
    }  // end loop over entities
  }    // end loop over co-dimensions
}
}  // namespace lf::mesh::utils
