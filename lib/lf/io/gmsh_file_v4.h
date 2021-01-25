/**
 * @file
 * @brief Functionality to read a GmshFile Version 4 into a memory data
 * structure and continue processing it later on.
 * @author Raffael Casagrande
 * @date   2019-08-19 06:05:01
 * @copyright MIT License
 */

#ifndef __0eae50bb4868430ebde754b70d7f7f47
#define __0eae50bb4868430ebde754b70d7f7f47

#include <lf/mesh/mesh.h>

#include <optional>

namespace lf::io {

/// A representation of a .msh file (V4) in a c++ data structure.
struct GMshFileV4 {
  /// The version of GMSH of the msh file, equals usually 4.1
  std::string version_number = "4.1";
  /// Is it a binary file?
  bool is_binary = false;
  /// how many bytes is a size_t?
  int size_t_size = 8;
  /**
   * \brief Represents a physical name as defined in gmsh.
   * In GMSH a Physical name is created through one of the commands
   * `Physical Point`, `Physical Line`, `Physical Surface` or `Physical Volume`
   * and represents a collection of points, lines, surfaces or volumes.
   *
   *  Every Physical Entity has:
   *  - A dimension (of the entities that it contains)
   *  - A (integer) tag which must be unique per dimension, i.e. a point and a
   *    line can have the same number
   *  - A name (string)
   *
   **/
  struct PhysicalName {
    /// Physical dimension of this physical name (1 for lines, 2 for surfaces
    /// etc.)
    int dimension = 0;
    /**
     *  \brief The identification number of the Physical name
     *  This number is assigned to a physical entity as the first/second
     * argument (if there is a name) argument of a `Physical Point`,
     * `Physical Line`, `Physical Surface` or `Physical Volume` command.
     */
    int physical_tag = 0;
    /// The name of this Physical Entity (provided
    std::string name;
  };

  /**
   * \brief A list of all Physical entities that have a name.
   * \note If the user creates a physical entity in GMSH without a name then it
   * is not listed in this collection.
   */
  std::vector<PhysicalName> physical_names;

  /**
   * @brief An entity of dimension 0. Usually an entity corresponds to a
   * geometric construct as it is defined in the original *.geo file. Every
   * entity consists of one or more elements and nodes. And every entity can
   * belong to some physical tags.
   */
  struct PointEntity {
    /// A unique number that identifies this entity (there may be another entity
    /// with the same number but other dimension).
    int tag = -1;
    /// the coordinates of the point
    Eigen::Vector3d coord;
    /// The physical tags to which this entity belongs.
    std::vector<int> physical_tags;
  };

  /**
   * @brief A higher-dimensional entity (which is not a point) such as a curve,
   * surface or volume. It usually corresponds to an entity as it is created in
   * a *.geo file. Every such entity consists of nodes and elements and is
   * bounded by lower-dimensional entities.
   *
   * @note A Gmsh entity is not the same as a LehrFEM++ entity! A LehrFEM++
   * entity is a concrete mesh element such as a triangle or a quadrilateral. A
   * Gmsh entity is a more abstract concept such as a circle or a surface. Every
   * gmsh entity consists of one or more mesh elements (e.g. triangles).
   */
  struct Entity {
    /// unique number that identifies this entity (entities of different
    /// dimensions can have the same number!)
    int tag = -1;
    /**
     * @brief `min_coord` and `max_coord` define a bounding box within which the
     * entity lies.
     */
    Eigen::Vector3d min_coord;
    /// @copydoc min_coord
    Eigen::Vector3d max_coord;

    /**
     * @brief The physical tags that identify the physical names to which this
     * entity belongs
     */
    std::vector<int> physical_tags;
    /**
     * @brief The bounding entities (one dimension smaller) that bound this
     * entity.
     */
    std::vector<int> bounding_entities;
  };

  /**
   * @brief the (Gmsh-) entities present in this mesh. An entity typically
   * corresponds to a geometrical object specified during the creation of the
   * mesh, e.g. a circle or a volume enclosed by surfaces.
   *
   * `std::get<0>(entities)`: Point entities
   * `std::get<1>(entities)`: Curve entities
   * `std::get<2>(entities)`: Surface entities
   * `std::get<3>(entities)`: volume entities
   */
  std::tuple<std::vector<PointEntity>, std::vector<Entity>, std::vector<Entity>,
             std::vector<Entity>>
      entities;

  struct GhostEntity {
    /**
     * @brief uniquely identifies this ghost entity
     */
    int tag;

    /**
     * @brief the partition for which the elements in this ghost entity are
     * ghost entities.
     */
    int partition;
  };

  /**
   * @brief When a mesh is partitioned, the entities are also partitioned.
   *        This is a node entity which is typically used to describe other
   *        partitioned entities (as bounding entities)
   */
  struct PartitionedPointEntity {
    /**
     * @brief unique number that identifies this partitioned point.
     */
    int tag = -1;

    /**
     * @brief The dimension of the parent entity to which this partitioned node
     * belongs
     */
    int parent_dim = -1;

    /**
     * @brief The tag of the parent entity to which this partitioned entity
     * belongs
     */
    int parent_tag = -1;

    /**
     * @brief to which partitions does this entity belong?
     */
    std::vector<int> partitions;

    /**
     * @brief the coordinates of this point
     */
    Eigen::Vector3d coord;

    /**
     * @brief the physical tags to which this partitioned node entity belongs
     */
    std::vector<int> physical_tags;
  };

  /**
   * @brief A higher dimensional partitioned entity
   */
  struct PartitionedEntity {
    /**
     * @brief unique number that identifies this partitioned entity.
     */
    int tag = -1;

    /**
     * @brief The dimension of the parent entity to which this partitioned
     * entity belongs
     */
    int parent_dim = -1;

    /**
     * @brief The tag of the parent entity to which this partitioned entity
     * belongs
     */
    int parent_tag = -1;

    /**
     * @brief to which partitions does this entity belong?
     */
    std::vector<int> partitions;

    /**
     * @brief the `min_coord` and `max_coord` define the (axis-aligned) bounding
     * box within which this partitioned entity lies.
     */
    Eigen::Vector3d min_coord;

    /// @copydoc min_coord
    Eigen::Vector3d max_coord;

    /**
     * @brief the physical tags to which this partitioned entity belongs
     */
    std::vector<int> physical_tags;

    /**
     * @brief The tags of the partitioned entities of one smaller dimension
     * that bound this partitioned entity.
     */
    std::vector<int> bounding_entities;
  };

  struct PartitionedEntities {
    /**
     * @brief total number of mesh partitions
     * @sa partitioned_entities
     */
    std::size_t num_partitions = 0;

    /**
     * @brief Ghost entities of this mesh
     *
     * Every ghost entity consists of a number of ghost elements (see
     * ghost_elements) and defines for which partition these elements are ghost
     * elements.
     */
    std::vector<GhostEntity> ghost_entities;

    /**
     * @brief a list of partitioned entities in this mesh.
     *
     * These entities are only present if the mesh creator has chosen to
     * partition the mesh. In such a case, the \ref entities entities that make
     * up the mesh are formally partitioned into smaller entities which are
     * listed here. The elements then refer to partitioned entities.
     *
     * - `std::get<0>(partitioned_entities)` are points
     * - `std::get<1>(partitioned_entities)` are curves
     * - `std::get<2>(partitioned_entities)` are surfaces
     * - `std::get<3>(partitioned_entities)` are volumes
     */
    std::tuple<std::vector<PartitionedPointEntity>,
               std::vector<PartitionedEntity>, std::vector<PartitionedEntity>,
               std::vector<PartitionedEntity>>
        partitioned_entities;
  };

  /**
   * @brief Information about partitioned entities (in case the mesh has been
   * partitioned)
   */
  PartitionedEntities partitioned_entities;

  struct NodeBlock {
    /**
     * @brief Dimension of the (partitioned) entity to which the nodes in this
     * block belong
     */
    int entity_dim = -1;

    /**
     * @brief The tag of the (partitioned) entity to which the nodes in this
     * block belong
     */
    int entity_tag = -1;

    /**
     * @brief Are the nodes in this block parametric?
     */
    bool parametric = false;

    /**
     * @brief a list of nodes in this block
     *
     * - `std::get<0>(nodes[i])` contains the tag of the i-th node
     * - `std::get<1>(nodes[i])` contains the coordinates of the i-th node
     *
     * @note Currently we don't support the parametric coordinates because it's
     * not clear from the documentation how they work...
     */
    std::vector<std::pair<std::size_t, Eigen::Vector3d>> nodes;
  };

  struct Nodes {
    /// Total number of nodes in the mesh
    std::size_t num_nodes = -1;

    /// Smallest node tag that exists
    std::size_t min_node_tag = 0;

    /// biggest node tag that exists
    std::size_t max_node_tag = 0;

    /**
     * @brief The nodes that make up this mesh organized in blocks.
     *
     */
    std::vector<NodeBlock> node_blocks;
  };

  Nodes nodes;

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
    TRIA6 = 9,     //!< 6-node second order triangle (3 nodes associated with
                   //!< the vertices and 3 with the edges).
    QUAD9 = 10,    //!< 9-node second order quadrangle (4 nodes associated with
                   //!< the vertices, 4 with the edges and 1 with the face).
    TET10 = 11,    //!< 10-node second order tetrahedron (4 nodes associated
                   //!< with the vertices and 6 with the edges).
    HEX27 = 12,    //!< 27-node second order hexahedron (8 nodes associated
                   //!< with the vertices, 12 with the edges, 6 with the faces
                   //!< and 1 with the volume).
    PRISM18 = 13,  //!< 18-node second order prism (6 nodes associated with
                   //!< the vertices, 9 with the edges and 3 with the
                   //!< quadrangular faces).
    PYRAMID14 = 14,  //!< 14-node second order pyramid (5 nodes associated
                     //!< with the vertices, 8 with the edges and 1 with the
                     //!< quadrangular face).
    POINT = 15,      //!< 1-node point
    QUAD8 = 16,      //!< 8-node second order quadrangle (4 nodes associated
                     //!< with the vertices and 4 with the edges).
    HEX20 = 17,      //!< 20-node second order hexahedron (8 nodes associated
                     //!< with the vertices and 12 with the edges).
    PRISM15 = 18,    //!< 15-node second order prism (6 nodes associated with
                     //!< the vertices and 9 with the edges).
    PYRAMID13 = 19,  //!< 13-node second order pyramid (5 nodes associated
                     //!< with the vertices and 8 with the edges).
    TRIA9 = 20,      //!< 9-node third order incomplete triangle (3 nodes
                     //!< associated with the vertices, 6 with the edges)
    TRIA10 = 21,     //!< 10-node third order triangle (3 nodes associated with
                     //!< the vertices, 6 with the edges, 1 with the face)
    TRIA12 = 22,     //!< 12-node fourth order incomplete triangle (3 nodes
                     //!< associated with the vertices, 9 with the edges)
    TRIA15 = 23,     //!< 15-node fourth order triangle (3 nodes associated
                     //!< with the vertices, 9 with the edges, 3 with the face)
    TRIA15_5 = 24,   //!< 15-node fifth order incomplete triangle (3 nodes
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
    TET35 = 30,   //!< 35-node fourth order tetrahedron (4 nodes associated
                  //!< with the vertices, 18 with the edges, 12 with the
                  //!< faces, 1 in the volume)
    TET56 = 31,   //!< 56-node fifth order tetrahedron (4 nodes associated
                  //!< with the vertices, 24 with the edges, 24 with the
                  //!< faces, 4 in the volume)
    HEX64 = 92,   //!< 64-node third order hexahedron (8 nodes associated
                  //!< with the vertices, 24 with the edges, 24 with the
                  //!< faces, 8 in the volume)
    HEX125 = 93   //!< 125-node fourth order hexahedron (8 nodes associated
                  //!< with the vertices, 36 with the edges, 54 with the
                  //!< faces, 27 in the volume)
  };

  /// Contains a list of all element types that are possible.
  static const std::vector<ElementType> AllElementTypes;

  /**
   * @brief Represents number of mesh elements (such as
   * triangles/quadrilaterals/points/lines) that share the same properties
   * @attention In comparison to LehrFEM++, in GMSH every "entity" is referred
   * to as an "element"
   */
  struct ElementBlock {
    /// Dimension of the elements in this block
    int dimension = -1;

    /**
     * @brief Tag of the entity to which elements in this block belong. If the
     * mesh is not partitioned, this references a `entities.tag`. Otherwise it
     * references a partitioned entity `partitioned_entities.tag`.
     */
    int entity_tag = -1;

    /// @brief type of the elements in this block
    ElementType element_type = ElementType::EDGE2;

    /**
     * @brief The elements in this block
     *
     * - `std::get<0>(elements[i])` is the unique tag of this element
     * - `std::get<1>(elements[i])` lists the node tags that make up this
     * element.
     */
    std::vector<std::pair<std::size_t, std::vector<std::size_t>>> elements;
  };

  struct Elements {
    /// Total number of elements in the mesh
    std::size_t num_elements = 0;

    /// minimum element tag
    std::size_t min_element_tag = 0;

    /// maximum element tag
    std::size_t max_element_tag = 0;

    /**
     * @brief A list of all Elements (Points,Lines,Surfaces or Volumes) present
     * in the *.msh file organized in blocks.
     *
     * @note  GMSH writes out all surface
     * elements if the mesh is 2D and all volume elements if the mesh is 3d.
     * Points, Lines (and Surfaces for 3D mesh) are not written to the file in
     * general. They only appear if they are part of some physical entity.
     */
    std::vector<ElementBlock> element_blocks;
  };

  /**
   * @brief Information about all mesh elements in this file.
   */
  Elements elements;

  /**
   * @brief Describes how 2 elementary entities are identified with each other
   * to represent periodic boundaries.
   *
   * It contains the dimension of the elementary entities, the elementary entity
   * numbers on the slave and master side and a list of nodes on the
   * slave/master side that should be identified with each other and that make
   * up the slave/master elementary entity.
   */
  struct PeriodicLink {
    /// Dimension of the elementary entities/elements that are coupled to each
    /// other.
    int dimension = 0;
    /// The (partitioned) entity number of the slave side
    int entity_tag_slave = 0;
    /// The (partitioned) entity number of the master side
    int entity_tag_master = 0;

    /// The transformation matrix to map the slave side to the master side (in
    /// homogeneous coordinates).
    std::optional<Eigen::Matrix4d> affine_transform;

    /**
     * \brief A List of nodes that pairs nodes on the slave side with nodes on
     *the master side.
     *
     * - `nodeMapping[i].first` is the tag of a node on the \e
     *slave side, which belongs to the (partitioned) entity `entity_tag_slave`.
     * - `nodeMapping[i].first` is the tag of a node on the \e
     *slave side, which belongs to the (partitioned) entity `entity_tag_master`.
     **/
    std::vector<std::pair<std::size_t, std::size_t>> node_mapping;
  };

  /// A List of Periodic definitions identifying elementary entities on the
  /// boundary with each other.
  std::vector<PeriodicLink> periodic_links;

  struct GhostElement {
    /// @brief uniquely identifies a mesh element (@sa element_blocks)
    std::size_t element_tag = 0;

    /// @brief defines the partition to which this element belongs (mainly)
    int partition_tag = -1;

    /// @brief references a ghost_entities by tag, if this mesh element is a
    /// ghost entity for that mesh partition.
    std::vector<int> ghost_partition_tags;
  };

  /**
   * @brief A list of ghost elements in the mesh.
   */
  std::vector<GhostElement> ghost_elements;
};

/**
 * @brief Get the Reference Element type of a GmshElement
 * @param et The GmshFile element type
 */
base::RefEl RefElOf(GMshFileV4::ElementType et);

/// Dimension of the GmshElement type
int DimOf(GMshFileV4::ElementType et);

/// Number of nodes that this element type has
int NumNodes(GMshFileV4::ElementType et);

/// Output the element type onto a stream:
std::ostream& operator<<(std::ostream& stream, GMshFileV4::ElementType et);

/**
 * @brief Read a GmshFile with format 4 and return it as an in-memory struct
 *
 * @param begin beginning of the file to parse (without header)
 * @param end end of the file to parse
 * @param version the exact version of the file, should be >=4
 * @param is_binary Is the file in binary format? (from header)
 * @param size_t_size sizeof(std::size_t) (from header)
 * @param one Representation of the number one if `is_binary==true` (from
 * header)
 * @param filename The name of the file that is being parsed (for better
 * diagnostics)
 */
GMshFileV4 ReadGmshFileV4(std::string::const_iterator begin,
                          std::string::const_iterator end,
                          const std::string& version, bool is_binary,
                          int size_t_size, int one,
                          const std::string& filename);

}  // namespace lf::io

#endif  // __0eae50bb4868430ebde754b70d7f7f47
