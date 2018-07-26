/**
 * @file
 * @brief Implementation of PrintInfo functions
 * @author Raffael Casagrande
 * @date   2018-07-01 01:33:21
 * @copyright MIT License
 */

#include "print_info.h"
#include "lf/geometry/geometry.h"

namespace lf::mesh::utils {

void PrintInfo(const Mesh &mesh, std::ostream &o) {
  using dim_t = Mesh::dim_t;
  using size_type = Mesh::size_type;

  const dim_t dim_mesh = mesh.DimMesh();
  const dim_t dim_world = mesh.DimWorld();
  o << "Mesh of dimension " << dim_mesh << ", ambient dimension " << dim_world
    << std::endl;

  // Loop over codimensions
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


// Print function for Entity object
void PrintInfo(const lf::mesh::Entity& e, std::ostream& stream){

    lf::base::RefEl e_ref_el = e.RefEl();
    int dim_ref_el = e_ref_el.Dimension();
    stream << "Derived type of entity: " << typeid(e).name() << std::endl;
    stream << "Type of reference element: " << e_ref_el << std::endl;


    if(Entity::output_ctrl_ == 0){ // > 0
        int co_dim_entity = e.Codim();
        stream << "Codimension of entity w.r.t. Mesh.dimMesh(): " << co_dim_entity << std::endl;

        // Geometry of entity
        const geometry::Geometry *e_geo_ptr = e.Geometry();
        const Eigen::MatrixXd &ref_el_corners(e_ref_el.NodeCoords()); // From PrintInfo(Mesh)
        //int dim_global = e_geo_ptr->DimGlobal();
        //stream << std::endl << e_geo_ptr->Global(ref_el_corners) << std::endl; // Causes segmentation fault


        for(int co_dim = dim_ref_el; co_dim > 0; co_dim--){
            int num_sub_ent = e_ref_el.NumSubEntities(co_dim);
            stream << "Codimension " << co_dim << " has " << num_sub_ent << " entities" << std::endl;

            if(Entity::output_ctrl_ == 0){ // > 10
                // Geometry of subentities
                for (const Entity &sub_ent : e.SubEntities(0)){ // (co_dim_entities)
                    lf::base::RefEl sub_ent_refel = sub_ent.RefEl();
                    stream << "Subentity, type: " << sub_ent_refel << std::endl;
                    const geometry::Geometry *sub_e_geo_ptr = sub_ent.Geometry();
                }
            }
        }

    } // if

} // PrintInfo


/*
// ?? void PrintInfo(const lf::mesh::Entity& e, const lf::mesh::Mesh& m, std::ostream &stream);
*/

}  // namespace lf::mesh::utils
