/**
 * @file
 * @brief Functionality to read a Gmsh file Version 2.0 into a in-memory data
 * structure for later processing.
 * @author Raffael Casagrande
 * @date   2019-08-19 03:23:45
 * @copyright MIT License
 */

#ifndef __1fd687af60b74bf6aad50d509ecbc4da
#define __1fd687af60b74bf6aad50d509ecbc4da

#include <lf/mesh/mesh.h>

namespace lf::io {

/// A representation of a .msh file (V2) in a c++ data structure.
struct GMshFileV2 {
  using size_type = mesh::Mesh::size_type;
  /// The version of GMSH of the msh file, equals usually 2.2
  std::string VersionNumber;
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
std::ostream& operator<<(std::ostream& stream, const GMshFileV2& mf);

/// Output the element type onto the console:
std::ostream& operator<<(std::ostream& stream, GMshFileV2::ElementType et);

/// Number of nodes that this element type has
unsigned int NumNodes(GMshFileV2::ElementType et);

/// Dimension of the GmshElement type
int DimOf(GMshFileV2::ElementType et);

/// Reference element type of a GmshElementType
base::RefEl RefElOf(GMshFileV2::ElementType et);

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
GMshFileV2 readGmshFileV2(std::string::const_iterator begin,
                          std::string::const_iterator end,
                          const std::string& version, bool is_binary,
                          int size_t_size, int one,
                          const std::string& filename);

}  // namespace lf::io

#endif  // __1fd687af60b74bf6aad50d509ecbc4da
