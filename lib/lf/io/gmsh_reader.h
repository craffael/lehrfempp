/***************************************************************************
 * LehrFEM++ - A simple C++ finite element library for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Declares the class GmshReader
 * @author Raffael Casagrande
 * @date   2018-06-29 04:51:59
 * @copyright MIT License
 */

#ifndef INCG7fedf7cf1a0246a98b2bf431cfa34da2
#define INCG7fedf7cf1a0246a98b2bf431cfa34da2
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <map>
#include <string>
#include <variant>
#include <vector>

#include "gmsh_file_v2.h"
#include "gmsh_file_v4.h"

namespace lf::io {

/**
 * @brief Reads a [Gmsh](http://gmsh.info/) `*.msh` file into a
 * mesh::MeshFactory and provides a link between mesh::Entity objects and
 * the gmsh's physical entities.
 *
 * In order to import the `*.msh` file successfully make sure that:
 * - Save the mesh in Gmsh's proprietary `*.msh` file format
 *   (`Version 2 ASCII`, `Version 2 Binary`, `Version 4 ASCII` or `Version 4
 *   Binary`)
 * - If you have specified physical entities in Gmsh and you want to export
 *   them, make sure you don't tick `Save all (ignore physical groups)` and
 *   make sure that every surface (2d) / volume (3d) belongs to at least one
 *   physical entity.
 *   (Otherwise the corresponding mesh elements are not exported by Gmsh)
 * - If you didn't specify any physical entities in Gmsh, you should tick
 *   `Save all (ignore physical groups)`.
 * - The GmshReader doesn't support the options "Save parametric coordinates"
 *   and "Save one file per partition".
 *
 * The import of meshes from Gmsh-generated mesh files is also discussed in
 * [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{rem:betlgmsh}.
 *
 * #### Sample usage:
 * @snippet gmsh_reader.cc usage
 *
 * #### Note about mesh completeness:
 * In principle, Gmsh can export meshes that are not "complete", i.e. meshes
 * containing entities with codim>0 that are not sub-entities of cell-entities.
 * Also, if Gmsh exports a higher-order mesh, it has to introduce
 * "auxilliary nodes" in order to define the higher-order mesh elements (these
 * auxilliary nodes are normally placed on edges, faces or in the interior of
 * cells).
 *
 * Since incomplete meshes are usually not useful in a FEM context, GmshReader
 * doesn't insert into the mesh auxilliary nodes or nodes not belonging to at
 * least one cell. Auxilliary nodes are however used in the construction of
 * higher order geometry objects.
 *
 */
class GmshReader {
 public:
  using size_type = mesh::Mesh::size_type;
  using dim_t = base::RefEl::dim_t;
  using GmshFileVariant = std::variant<GMshFileV2, GMshFileV4>;

  /**
   * @brief Get the mesh that was read by this reader.
   */
  [[nodiscard]] std::shared_ptr<mesh::Mesh> mesh() { return mesh_; }

  /**
   * @brief Get the mesh that was read by this reader.
   */
  [[nodiscard]] std::shared_ptr<const mesh::Mesh> mesh() const { return mesh_; }

  /**
   * @brief maps the name of a physical entity to the physical entity number
   *        of gmsh.
   * @param name Name of the physical entity.
   * @param codim Optionally the codimension of the physical entity (w.r.t
   * dimMesh), only needed if there are physical entities with the same name but
   * different dimensions.
   * @return  Number of the physical entity (as shown in Gmsh)
   * @note If you have defined a physical entity in GMSH without
   * giving it a name, you cannot use this function.
   * @sa PhysicalEntityNr2Name() and the code snippet in the documentation of
   * @ref GmshReader.
   */
  [[nodiscard]] size_type PhysicalEntityName2Nr(const std::string& name,
                                                dim_t codim = -1) const;

  /**
   * @brief Gives the name of a physical entity (inverse of
   * PhysicalEntityName2Nr())
   * @param  number   The number of the physical entity (as obtained from
   * PhysicalEntityName2Nr() )
   * @param codim     Specify additionally the codimension of the physical
   * entity. If not specified and there is exactly one physical entity with the
   * given number, that entity is returned. If there are multiple entities with
   * the same number and `codim=-1`, an error is thrown.
   * @return The name of the physical entity with number `number`
   * @sa PhysicalEntityNr2Name
   */
  [[nodiscard]] std::string PhysicalEntityNr2Name(size_type number,
                                                  dim_t codim = -1) const;

  /**
   * @brief Retrieve a list of all (Gmsh) physical entities of the given codim.
   * @param codim The codimension
   * @return A list of physical entities (number, name)
   *
   * See the code snippet in the documentation of @ref GmshReader.
   */
  [[nodiscard]] std::vector<std::pair<size_type, std::string>> PhysicalEntities(
      dim_t codim) const;

  /**
   * \brief Return the Physical Entity Numbers associated with this mesh entity,
   *        sorted ascending.
   * \param e  The entity of the grid.
   * \return   The Physical Entity Number that was assigned in GMSH.
   */
  [[nodiscard]] std::vector<size_type> PhysicalEntityNr(
      const mesh::Entity& e) const;

  /**
   * @brief Test whether the given entity belongs to a Gmsh physical entity.
   * @param e The entity that should be tested.
   * @param physical_entity_nr The number of the gmsh physical entity.
   * @return True if the entity `e` belongs to the physical entity.
   *
   * This method is related to the method PhysicalEntityNr(): It returns true,
   * if the vector returned by `PhysicalEntityNr(e)` contains
   * `physical_entity_nr`
   *
   * @sa PhysicalEntityNr() and the code snippet in the documentation of
   * @ref GmshReader.
   */
  [[nodiscard]] bool IsPhysicalEntity(const mesh::Entity& e,
                                      size_type physical_entity_nr) const;

  /**
   * @brief Create a new GmshReader from the given MshFile (advanced usage)
   * @param factory The mesh::MeshFactory that is used to construct the mesh.
   * @param msh_file A LehrFEM++ specific representation of a GmshFile.
   *
   * @sa MshFile
   */
  GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
             const GmshFileVariant& msh_file);

  /**
   * @brief Create a new GmshReader by reading from the specified file.
   * @param factory The mesh::MeshFactory that is used to construct the mesh.
   * @param filename The filename of the `.msh` file that is read.
   *
   * @note If the `factory.DimWorld() == 3`, there must be at least one
   *       3D mesh element in the *.msh file. Similarly, if
   *       `factory.DimWorld() == 2` there should be only 2D mesh elements
   *       in the *.msh file!
   * @note GmshReader supports ASCII and Binary `.msh` files.
   */
  GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
             const std::string& filename);

 private:
  /// The underlying grid created by the grid factory.
  std::shared_ptr<mesh::Mesh> mesh_;

  std::unique_ptr<mesh::MeshFactory> mesh_factory_;

  /// The PhysicalEntityNr of every entity (0 if not set):
  std::shared_ptr<mesh::utils::AllCodimMeshDataSet<std::vector<size_type>>>
      physical_nrs_;

  /// Map from physicalEntity name -> nr, codim
  std::multimap<std::string, std::pair<size_type, dim_t>> name_2_nr_;

  /// Map from physicalEntity nr -> name, codim
  std::multimap<size_type, std::pair<std::string, dim_t>> nr_2_name_;

  void InitGmshFile(const GMshFileV2& msh_file);
  void InitGmshFile(const GMshFileV4& msh_file);
};

std::variant<GMshFileV2, GMshFileV4> ReadGmshFile(const std::string& filename);

}  // namespace lf::io

#endif  // INCG7fedf7cf1a0246a98b2bf431cfa34da2
