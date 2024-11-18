/**
 * @file
 * @brief Implementation of PrintInfo functions
 * @author Raffael Casagrande
 * @date   2018-07-01 01:33:21
 * @copyright MIT License
 */

#include "print_info.h"

#include <typeinfo>

#include "lf/geometry/geometry.h"

namespace lf::mesh::utils {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void PrintInfo(std::ostream &o, const lf::mesh::Mesh &mesh, int ctrl) {
  using dim_t = lf::base::dim_t;
  using size_type = lf::base::size_type;

  const dim_t dim_mesh = mesh.DimMesh();
  const dim_t dim_world = mesh.DimWorld();
  o << "Mesh of dimension " << static_cast<int>(dim_mesh)
    << ", ambient dimension " << static_cast<int>(dim_world) << '\n';

  if (ctrl > 10) {
    // Loop over codimensions

    for (int co_dim = base::narrow<int>(dim_mesh); co_dim >= 0; co_dim--) {
      const size_type no_ent = mesh.NumEntities(co_dim);
      o << "Co-dimension " << co_dim << ": " << no_ent << " entities" << '\n';

      // Loop over entities
      for (const Entity *e : mesh.Entities(co_dim)) {
        const size_type e_idx = mesh.Index(*e);
        const dim_t e_codim = e->Codim();
        const geometry::Geometry *e_geo_ptr = e->Geometry();
        const lf::base::RefEl e_refel = e->RefEl();

        LF_VERIFY_MSG(e_geo_ptr,
                      co_dim << "-entity " << e_idx << ": missing geometry");
        LF_VERIFY_MSG(e_geo_ptr->DimLocal() == dim_mesh - co_dim,
                      co_dim << "-entity " << e_idx << ": wrong dimension");
        LF_VERIFY_MSG(e_geo_ptr->RefEl() == e_refel,
                      co_dim << "-entity " << e_idx << ": refEl mismatch");
        LF_VERIFY_MSG(e_codim == co_dim, co_dim << "-entity " << e_idx
                                                << " co-dimension mismatch");
        const Eigen::MatrixXd &ref_el_corners(e_refel.NodeCoords());
        o << "entity " << e_idx << " (" << e_refel << "): ";

        // Loop over local co-dimensions
        if (ctrl > 90) {
          for (int l = 1; l <= base::narrow<int>(dim_mesh) - co_dim; l++) {
            o << "rel codim-" << l << " subent: [";
            // Fetch subentities of co-dimension l
            auto sub_ent_range = e->SubEntities(l);
            for (const Entity *sub_ent : sub_ent_range) {
              o << mesh.Index(*sub_ent) << ' ';
            }
            o << "]";
          }
        }
        if (ctrl > 50) {
          o << '\n' << e_geo_ptr->Global(ref_el_corners);
        }
        o << '\n';
      }  // end loop over entities
    }  // end loop over co-dimensions
  }  // end printinfo_ctrl > 10
}  // end function PrintInfo

// Print function for Entity object
void PrintInfo(std::ostream &stream, const lf::mesh::Entity &e,
               int output_ctrl) {
  // Topological type of entity
  const lf::base::RefEl e_ref_el = e.RefEl();
  const base::dim_t dim_ref_el = e_ref_el.Dimension();
  // Geometry of entity and coordinates
  const geometry::Geometry *e_geo_ptr = e.Geometry();
  LF_ASSERT_MSG(e_geo_ptr != nullptr, "Missing geometry information!");

  stream << "Entity " << e_ref_el << "/" << typeid(*e_geo_ptr).name() << '\n';

  if (output_ctrl > 10) {
    stream << "Dimension: " << dim_ref_el << '\n';

    const Eigen::MatrixXd &ref_el_corners(e_ref_el.NodeCoords());

    // Loop over codimensions
    for (base::dim_t co_dim = dim_ref_el; co_dim > 0; co_dim--) {
      const base::size_type num_sub_ent = e_ref_el.NumSubEntities(co_dim);
      stream << '\n'
             << "Codimension " << co_dim << ": " << num_sub_ent
             << " sub-entities" << '\n';

      if (output_ctrl > 50) {
        int sub_ent_num = 0;
        // Loop over subentities
        for (const Entity *sub_ent : e.SubEntities(co_dim)) {
          const lf::base::RefEl sub_ent_refel = sub_ent->RefEl();
          stream << "* Subentity " << sub_ent_num << " (" << sub_ent_refel
                 << ")" << '\n';

          if (output_ctrl > 90) {
            // Print coordinates
            stream << e_geo_ptr->Global(ref_el_corners).col(sub_ent_num)
                   << '\n';
          }

          sub_ent_num += 1;

        }  // loop sub-ent
      }  // if output ctrl
    }  // loop codim

    if (e_ref_el == lf::base::RefEl::kPoint() && output_ctrl > 90) {
      stream << e_geo_ptr->Global(ref_el_corners) << '\n';
    }

    stream << "-----------------------" << '\n';

  }  // if
}  // PrintInfo

}  // end namespace lf::mesh::utils
