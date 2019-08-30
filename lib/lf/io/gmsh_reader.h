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

#ifndef __7fedf7cf1a0246a98b2bf431cfa34da2
#define __7fedf7cf1a0246a98b2bf431cfa34da2
#include <lf/mesh/mesh.h>
#include <map>
#include <string>
#include <vector>
#include "lf/mesh/utils/utils.h"

namespace lf::io {
/// A representation of a .msh file in a c++ data structure.
struct MshFile {
  using size_type = mesh::Mesh::size_type;
  /// The version of GMSH of the msh file, equals usually 2.2
  double VersionNumber = 0;
  /// Is it a binary file?
  bool IsBinary = false;
  /// how many bytes is a double?
  int DoubleSize = 64;

  /**
   * \brief Represents a physical entity as defined in gmsh.
   * In GMSH a Physical entity is created through one of the commands `Physical
   *Point`, `Physical Line`, `Physical Surface` or `Physical Volume` and
   *represents a collection of points, lines, surfaces or volumes.
   *
   *  Every Physical Entity has:
   *  - A dimension (of the entities that it contains)
   *  - A Number which must be unique per dimension, i.e. a point and a line can
   *have the same number
   *  - A name (string)
   *
   **/
  struct PhysicalEntity {
    /// Physical dimension of this physical name (1 for lines, 2 for surfaces
    /// etc.)
    int Dimension = 0;
    /// The identification number of the physical entity (is specified in GMSH
    /// with the command `Physical Point`, `Physical Line`, Surface etc.)
    /**
     *  \brief The identification number of the Physical Entity
     *  This number is assigned to a physical entity as the first/second (if
     * there is a name) argument of a `Physical Point`, `Physical Line`,
     * `Physical Surface` or `Physical Volume` command.
     */
    int Number = 0;
    /// The name of this Physical Entity (provided
    std::string Name;
  };

  /**
   * \brief A list of all Physical entities that have a name.
   * \note If the user creates a physical entity in GMSH without a name then it
   * is not listed in this collection.
   */
  std::vector<PhysicalEntity> PhysicalEntities;

  /**
   * \brief The nodes that make up this mesh.
   *
   * - `node[i].first` contains the node number assigned by GMSH. They are not
   * necessarily dense or ordered in sequence.
   * - `node[i].second` contains the 3D coordinates of this point, for 2D mesh
   * the z-component is always zero.
   */
  std::vector<std::pair<size_type, Eigen::Vector3d>> Nodes;

  /// All possible element types (see GMSH documentation)
  enum class ElementType : int {
    EDGE2 = 1,     //!< 2-node line
    TRIA3 = 2,     //!< 3-node triangle
    QUAD4 = 3,     //!< 4-node quadrangle
    TET4 = 4,      //!< 4-node tetrahedron
    HEX8 = 5,      //!< 8-node hexahedron
    PRISM6 = 6,    //!< 6-node prism
    PYRAMID5 = 7,  //!< 5-node pyramid
    EDGE3 = 8,     //!< 3-node second order line (2 nodes associated with the
                   //!< vertices and 1 with the edge)
    TRIA6 = 9,   //!< 6-node second order triangle (3 nodes associated with the
                 //!< vertices and 3 with the edges).
    QUAD9 = 10,  //!< 9-node second order quadrangle (4 nodes associated with
                 //!< the vertices, 4 with the edges and 1 with the face).
    TET10 = 11,  //!< 10-node second order tetrahedron (4 nodes associated with
                 //!< the vertices and 6 with the edges).
    HEX27 = 12,  //!< 27-node second order hexahedron (8 nodes associated with
                 //!< the vertices, 12 with the edges, 6 with the faces and 1
                 //!< with the volume).
    PRISM18 =
        13,  //!< 18-node second order prism (6 nodes associated with the
             //!< vertices, 9 with the edges and 3 with the quadrangular faces).
    PYRAMID14 =
        14,  //!< 14-node second order pyramid (5 nodes associated with the
             //!< vertices, 8 with the edges and 1 with the quadrangular face).
    POINT = 15,    //!< 1-node point
    QUAD8 = 16,    //!< 8-node second order quadrangle (4 nodes associated with
                   //!< the vertices and 4 with the edges).
    HEX20 = 17,    //!< 20-node second order hexahedron (8 nodes associated with
                   //!< the vertices and 12 with the edges).
    PRISM15 = 18,  //!< 15-node second order prism (6 nodes associated with the
                   //!< vertices and 9 with the edges).
    PYRAMID13 = 19,  //!< 13-node second order pyramid (5 nodes associated with
                     //!< the vertices and 8 with the edges).
    TRIA9 = 20,  //!< 9-node third order incomplete triangle (3 nodes associated
                 //!< with the vertices, 6 with the edges)
    TRIA10 = 21,  //!< 10-node third order triangle (3 nodes associated with the
                  //!< vertices, 6 with the edges, 1 with the face)
    TRIA12 = 22,  //!< 12-node fourth order incomplete triangle (3 nodes
                  //!< associated with the vertices, 9 with the edges)
    TRIA15 = 23,  //!< 15-node fourth order triangle (3 nodes associated with
                  //!< the vertices, 9 with the edges, 3 with the face)
    TRIA15_5 = 24,  //!< 15-node fifth order incomplete triangle (3 nodes
                    //!< associated with the vertices, 12 with the edges)
    TRIA21 = 25,  //!< 21-node fifth order complete triangle (3 nodes associated
                  //!< with the vertices, 12 with the edges, 6 with the face)
    EDGE4 = 26,   //!< 4-node third order edge (2 nodes associated with the
                  //!< vertices, 2 internal to the edge)
    EDGE5 = 27,   //!< 5-node fourth order edge (2 nodes associated with the
                  //!< vertices, 3 internal to the edge)
    EDGE6 = 28,   //!< 6-node fifth order edge (2 nodes associated with the
                  //!< vertices, 4 internal to the edge)
    TET20 = 29,   //!< 20-node third order tetrahedron (4 nodes associated with
                  //!< the vertices, 12 with the edges, 4 with the faces)
    TET35 = 30,   //!< 35-node fourth order tetrahedron (4 nodes associated with
                  //!< the vertices, 18 with the edges, 12 with the faces, 1 in
                  //!< the volume)
    TET56 = 31,   //!< 56-node fifth order tetrahedron (4 nodes associated with
                  //!< the vertices, 24 with the edges, 24 with the faces, 4 in
                  //!< the volume)
    HEX64 = 92,   //!< 64-node third order hexahedron (8 nodes associated with
                  //!< the vertices, 24 with the edges, 24 with the faces, 8 in
                  //!< the volume)
    HEX125 = 93   //!< 125-node fourth order hexahedron (8 nodes associated with
                  //!< the vertices, 36 with the edges, 54 with the faces, 27 in
                  //!< the volume)
  };

  /// Contains a list of all element types that are possible.
  static const std::vector<ElementType> AllElementTypes;

  /**
   * \brief Represents a mesh volume/surface/line/point
   * \attention In comparison to HyDi in GMSH every "entity" is referred to as
   * an "element"
   */
  struct Element {
    /// The number of this element in the mesh. They are not necessarily dense
    /// or ordered in sequence.
    size_type Number = 0;
    /// The element type
    ElementType Type = ElementType::POINT;
    /**
     * \brief The Number of the Physical Entity to which this element belongs
     * (this is the first tag written in the .msh file) \note  Two Physical
     * Entities with different dimensions can have the same number!
     */
    int PhysicalEntityNr = 0;
    /**
     * \brief The number of the elementary entity to which this element belongs
     * (second element tag in .msh file)
     *
     * An elementary entity is an entity that was specifically set by the user
     * to define the geometry that should be meshed. This is the number that was
     * given by the user to the point/line/surface/volume.
     */
    int ElementaryEntityNr = 0;

    /**
     * \brief The id's of the partition to which this element belongs.
     * \note Negative partition id's indicate ghost cells.
     * \note If the user did not partition the mesh, this is just an empty
     * vector.
     */
    std::vector<int> MeshPartitions;

    /**
     * \brief Contains the node numbers that make up this element (depends on
     * the element type to) \note The node number in this vector agrees with
     * `Nodes.first`.
     */
    std::vector<size_type> NodeNumbers;
  };

  /**
   * \brief A list of all Elements (Points,Lines,Surfaces or Volumes) present in
   * the *.msh file. \note  GMSH writes out all surface elements if the mesh is
   * 2D and all volume elements if the mesh is 3d. Points, Lines (and Surfaces
   * for 3D mesh) are not written to the file in general. They only appear if
   * they are part of some physical entity.
   */
  std::vector<Element> Elements;

  /**
   * \brief Describes how 2 elementary entities are identified with each to
   * represent periodic boundaries.
   *
   * It contains the dimension of the elementary entities, the elementary entity
   * numbers on the slave and master side and a list of nodes on the
   * slave/master side that should be identified with each other and that make
   * up the slave/master elementary entity.
   */
  struct PeriodicEntity {
    /// Dimension of the elementary entities that are coupled to each other.
    int Dimension = 0;
    /// The elementary entity number (\sa Element) on the slave side.
    int ElementarySlaveNr = 0;
    /// The elementary entity number (\sa Element) on the master side.
    int ElementaryMasterNr = 0;
    /**
     * \brief A List of nodes that pairs nodes on the slave side with nodes on
     *the master side.
     *
     * - `nodeMapping[i].first` is the number of a node (\sa Nodes) on the \e
     *slave side, which belongs to the elementary entity `ElementarySlaveNr`.
     * - `nodeMapping[i].second` is the number of a node (\sa Nodes) on the \e
     *master side, which belongs to the elementary entity `ElementaryMasterNr`
     **/
    std::vector<std::pair<size_type, size_type>> NodeMapping;
  };

  /// A List of Periodic definitions identifying elementary entities on the
  /// boundary with each other.
  std::vector<PeriodicEntity> Periodic;
};

/// For debugging purposes: Write the MshFile into a stream
std::ostream& operator<<(std::ostream& stream, const MshFile& mf);

/// Output the element type onto the console:
std::ostream& operator<<(std::ostream& stream, MshFile::ElementType et);

/// Number of nodes that this element type has
unsigned int NumNodes(MshFile::ElementType et);

/// Dimension of the GmshElement type
int DimOf(MshFile::ElementType et);

/**
 * \brief Read a *.msh file from disk and copy it's contents into the MshFile
 * Datastructure.
 *
 * So far the following sections of the .msh file are read:
 * - `$MeshFormat`
 * - `PhysicalNames`
 * - `$Nodes`
 * - `$Elements`
 * - `$Periodic`
 *
 * All other sections are ignored.
 *
 * \note We support the MshFile format 2.2 in binary or text form.
 * \note This routine is mainly used by the GmshReader class.
 */
MshFile readGmshFile(std::string path);

/**
 * @brief Reads a [Gmsh](http://gmsh.info/) `*.msh` file into a
 * mesh::MeshFactory and provides a link between mesh::Entity objects and
 * the gmsh's physical entities.
 *
 * In order to import the `*.msh` file successfully make sure that:
 * - Save the mesh in Gmsh's proprietary `*.msh` file format
 *   (`Version 2 ASCII` or `Version 2 Binary`)
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
 *
 *
 */
class GmshReader {
 public:
  using size_type = mesh::Mesh::size_type;
  using dim_t = base::RefEl::dim_t;

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
   * @sa PhysicalEntityNr2Name()
   */
  [[nodiscard]] size_type PhysicalEntityName2Nr(const std::string& name,
                                                dim_t codim = -1) const;

  /**
   * @brief Gives the name of a physical entity (inverse of
   * PhysicalEntityName2Nr())
   * @param  number   The number of the physical entity (as obtained from
   * PhysicalEntityName2Nr() )
   * @return The name of the physical entity with number `number`
   * @sa PhysicalEntityNr2Name
   */
  [[nodiscard]] std::string PhysicalEntityNr2Name(size_type number,
                                                  dim_t codim = -1) const;

  /**
   * @brief Retrieve a list of all (Gmsh) physical entities of the given codim.
   * @param codim The codimension
   * @return A list of physical entities (number, name)
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
   * @sa PhysicalEntityNr()
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
             const MshFile& msh_file);

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

  /// The PhysicalEntityNr of every node (0 if not set):
  std::shared_ptr<mesh::utils::AllCodimMeshDataSet<std::vector<size_type>>>
      physical_nrs_;

  /// Map from physicalEntity name -> nr, codim
  std::multimap<std::string, std::pair<size_type, dim_t>> name_2_nr_;

  /// Map from physicalEntity nr -> name, codim
  std::multimap<size_type, std::pair<std::string, dim_t>> nr_2_name_;
};

}  // namespace lf::io

#endif  // __7fedf7cf1a0246a98b2bf431cfa34da2
