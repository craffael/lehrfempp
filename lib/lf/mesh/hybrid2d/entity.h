#ifndef __7a3f1903d42141a3b1135e8e5ad72c1c
#define __7a3f1903d42141a3b1135e8e5ad72c1c

#include <lf/base/base.h>
#include "lf/mesh/entity.h"
#include "lf/mesh/mesh_interface.h"

namespace lf::mesh::hybrid2d {

class Mesh;

/**
 * @brief classes for topological entities in a 2D hybrid mesh
 * @{
 */

/**
 * @brief A node object for a 2D hybrid mesh
 * @tparam CODIM the co-dimension of the entity object \f$\in\{0,1,2\}\f$
 * 
 * @note Every `Entity` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Node : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Node() = default;

  /** @ brief Default and disabled constructors
   * @{ */
  Node(const Node&) = delete;
  Node(Node&&) noexcept = default;
  Node& operator=(const Node&) = delete;
  Node& operator=(Node&&) noexcept = default;
  /** @} */
  
 /**
  * @brief constructor, is called from MeshFactory
  * @param index index of the entity to be created; will usually be
  * retrieved via the `Index()` method of `Mesh`
  * @param geometry pointer to a geometry object providing the shape of the
  * entity
  *
  * @note Note that you need to create a suitable geometry object for the
  * entity before you can initialize the entity object itseld.
  */
  explicit Node(size_type index,
		std::unique_ptr<geometry::Geometry>&& geometry):
    index_(index),geometry_(std::move(geometry)) {
    LF_VERIFY_MSG(geometry->DimLocal() == 0,"Geometry must be that of a point");
    LF_VERIFY_MSG(geometry->RefEl() == base::RefEl::kPoint(),
		  "Geometry must fit point");
  }

  char Codim() const override { return 2; }

  base::RandomAccessRange<const mesh::Entity>
  SubEntities(char rel_codim) const override {
    LF_ASSERT_MSG(rel_codim == 0,"A point has only codim = 0 sub-entities");
    /* 
       return a range comprising only 'this'
    */
  }

  geometry::Geometry* Geometry() const override { return geometry_.get(); }

  base::RefEl RefEl() const override { return base::RefEl::kPoint(); }

  bool operator==(const mesh::Entity& rhs) const override { return this == &rhs; }

  virtual ~Node() override = default;

 private:
  size_type index_ = -1;                         // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_; // shape information
};

 /**
 * @brief An edge object for a 2D hybrid mesh
 * 
 * An topological edge object is define through two distinct references 
 * to node objects of the mesh. Their ordering reflects the intrinsic 
 * orientation of the mesh
 * @note Every `Edge` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Edge : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Edge() = default;

  /** @defgroup Default and disabled constructors
   * @{ */
  Edge(const Edge&) = delete;
  Edge(Edge&&) noexcept = default;
  Edge& operator=(const Edge&) = delete;
  Edge& operator=(Edge&&) noexcept = default;
  /** @} */
  
 /**
  * @brief constructor, is called from MeshFactory
  * @param index index of the entity to be created; will usually be
  * retrieved via the `Index()` method of `Mesh`
  * @param geometry pointer to a geometry object providing the shape of the edge
  * @param endpoint0 pointer to the first node
  * @param endpoint1 pointer to the second node
  *
  * @note Note that you need to create a suitable geometry object for the
  * entity before you can initialize the entity object itseld.
  */
  explicit Edge(size_type index,
		std::unique_ptr<geometry::Geometry>&& geometry,
		const Node *endpoint0,const Node *endpoint1):
    index_(index),geometry_(std::move(geometry)),endnodes_({endpoint0,endpoint1}) {
    LF_VERIFY_MSG((endpoint0 != nullptr) && (endpoint1 != nullptr),
		  "Invalid pointer to endnode of edge");
    LF_VERIFY_MSG(geometry->DimLocal() == 1,"Geometry must describe a curve");
    LF_VERIFY_MSG(geometry->RefEl() == base::RefEl::kSegment(),
		  "Edge geometry must fit a segment");
  }

  /** @brief an edge is an entity of co-dimension 1 */
  char Codim() const override { return 1; }

  /** @brief Access to all subentities selected by **relative** co-dimension 
   * @param rel_codim if 1 select endnodes, if 0 select edge itself
   * @return 
     - for rel_codim == 1: return 2-range covering endnodes
     - for rel_codim == 0: return the Edge entity itself
   */
  base::RandomAccessRange<const mesh::Entity>
    SubEntities(char rel_codim) const override;

  /** @defgroup Standard methods of an Entity object 
   * @sa mesh::Entity
   * @{
   */
  geometry::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return base::RefEl::kSegment(); }
  bool operator==(const mesh::Entity& rhs) const override { return this == &rhs; }
  /** @} */
  
  virtual ~Edge() override = default;
 private:
  size_type index_ = -1;                         // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_; // shape information
  std::array<const Node *,2> endnodes_;            // nodes connected by edge
};

 /**
 * @brief Describes a trilateral cell for a 2D hybrid mesh
 * 
 * A trilateral cell is defined by ordered lists of references to its nodes
 * and its edges; internal consistency is required
 * @note Every `Edge` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Trilateral : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Trilateral() = default;

  /** @defgroup Default and disabled constructors
   * @{ */
  Trilateral(const Trilateral&) = delete;
  Trilateral(Trilateral&&) noexcept = default;
  Trilateral& operator=(const Trilateral&) = delete;
  Trilateral& operator=(Trilateral&&) noexcept = default;
  /** @} */
  
 /**
  * @brief constructor, is called from MeshFactory
  * @param index index of the entity to be created; will usually be
  * retrieved via the `Index()` method of `Mesh`
  * @param geometry pointer to a geometry object providing the shape of the cell
  * @param corner0 pointer to first node
  * @param corner1 pointer to second node
  * @param corner2 pointer to third node
  * @param edge0 pointer to first edge
  * @param edge1 pointer to second edge
  * @param edge2 pointer to third edge
  *
  * The sub-entities have to be consistent according to the conventions
  * fixed for a reference element of type `kTria`.
  * @note Note that you need to create a suitable geometry object for the
  * entity before you can initialize the entity object itseld.
  */
  explicit Trilateral(size_type index,
		std::unique_ptr<geometry::Geometry>&& geometry,
		const Node *corner0,const Node *corner1,const Node *corner2,
		const Edge *edge0,const Edge *edge1,const Edge *edge2):
    index_(index),geometry_(std::move(geometry)),
    corners_({corner0,corner1,corner2}),edges_({edge0,edge1,edge2})
  {
    LF_VERIFY_MSG(corner0 != nullptr,"Invalid pointer to corner 0");
    LF_VERIFY_MSG(corner1 != nullptr,"Invalid pointer to corner 1");
    LF_VERIFY_MSG(corner2 != nullptr,"Invalid pointer to corner 2");
    LF_VERIFY_MSG(edge0 != nullptr,"Invalid pointer to edge 0");
    LF_VERIFY_MSG(edge1 != nullptr,"Invalid pointer to edge 1");
    LF_VERIFY_MSG(edge2 != nullptr,"Invalid pointer to edge 2");
    LF_VERIFY_MSG(geometry->DimLocal() == 2,"Geometry must describe a 2D cell");
    LF_VERIFY_MSG(geometry->RefEl() == base::RefEl::kTria(),
		  "Cell geometry must fit a triangle");
    /* 
       TODO: consistency check
    */
  }

  /** @brief an edge is an entity of co-dimension 1 */
  char Codim() const override { return 0; }

  /** @brief Access to all subentities selected by **relative** co-dimension 
   * @param rel_codim if 1 select edges, if 2 select nodes, if 0 select cell itself
   * @return 
     - for rel_codim == 1: return 3-range covering corners
     - for rel_codim == 2: return 3-range containing the edges
   */
  base::RandomAccessRange<const mesh::Entity>
    SubEntities(char rel_codim) const override;

  /** @defgroup Standard methods of an Entity object 
   * @sa mesh::Entity
   * @{
   */
  geometry::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return base::RefEl::kTria(); }
  bool operator==(const mesh::Entity& rhs) const override { return this == &rhs; }
  /** @} */
  
  virtual ~Trilateral() override = default;
 private:
  size_type index_ = -1;                         // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_; // shape information
  std::array<const Node *,3> corners_;           // nodes = corners of cell
  std::array<const Edge *,3> edges_;             // edges of the cells
};

 /**
 * @brief Describes a quadrilateral cell for a 2D hybrid mesh
 * 
 * A quadrilateral cell stores and ordered list of four nodes and
 * of four edges, which have to be compatible. 
 * @note Every `Edge` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Quadrilateral : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Quadrilateral() = default;

  /** @defgroup Default and disabled constructors
   * @{ */
  Quadrilateral(const Quadrilateral&) = delete;
  Quadrilateral(Quadrilateral&&) noexcept = default;
  Quadrilateral& operator=(const Quadrilateral&) = delete;
  Quadrilateral& operator=(Quadrilateral&&) noexcept = default;
  /** @} */
  
 /**
  * @brief constructor, is called from MeshFactory
  * @param index index of the entity to be created; will usually be
  * retrieved via the `Index()` method of `Mesh`
  * @param geometry pointer to a geometry object providing the shape of the cell
  * @param corner0 pointer to first node
  * @param corner1 pointer to second node
  * @param corner2 pointer to third node
  * @param edge0 pointer to first edge
  * @param edge1 pointer to second edge
  * @param edge2 pointer to third edge
  *
  * The sub-entities have to be consistent according to the conventions
  * fixed for a reference element of type `kTria`.
  * @note Note that you need to create a suitable geometry object for the
  * entity before you can initialize the entity object itseld.
  */
  explicit Quadrilateral(size_type index,
		std::unique_ptr<geometry::Geometry>&& geometry,
			 const Node *corner0,const Node *corner1,
			 const Node *corner2,const Node *corner3,
			 const Edge *edge0,const Edge *edge1,
			 const Edge *edge2,const Edge *edge3):
    index_(index),geometry_(std::move(geometry)),
    corners_({corner0,corner1,corner2,corner3}),edges_({edge0,edge1,edge2,edge3})
  {
    LF_VERIFY_MSG(corner0 != nullptr,"Invalid pointer to corner 0");
    LF_VERIFY_MSG(corner1 != nullptr,"Invalid pointer to corner 1");
    LF_VERIFY_MSG(corner2 != nullptr,"Invalid pointer to corner 2");
    LF_VERIFY_MSG(corner3 != nullptr,"Invalid pointer to corner 3");
    LF_VERIFY_MSG(edge0 != nullptr,"Invalid pointer to edge 0");
    LF_VERIFY_MSG(edge1 != nullptr,"Invalid pointer to edge 1");
    LF_VERIFY_MSG(edge2 != nullptr,"Invalid pointer to edge 2");
    LF_VERIFY_MSG(edge3 != nullptr,"Invalid pointer to edge 3");
    LF_VERIFY_MSG(geometry->DimLocal() == 2,"Geometry must describe a 2D cell");
    LF_VERIFY_MSG(geometry->RefEl() == base::RefEl::kQuad(),
		  "Cell geometry must fit a quad");
    /* 
       TODO: consistency check
    */
  }

  /** @brief an edge is an entity of co-dimension 1 */
  char Codim() const override { return 0; }

  /** @brief Access to all subentities selected by **relative** co-dimension 
   * @param rel_codim if 1 select edges, if 2 select nodes, if 0 select cell itself
   * @return 
     - for rel_codim == 1: return range with quadrilateral itself as only element
     - for rel_codim == 1: return 4-range covering corners
     - for rel_codim == 2: return 4-range containing the edges
   */
  base::RandomAccessRange<const mesh::Entity>
    SubEntities(char rel_codim) const override;

  /** @defgroup Standard methods of an Entity object 
   * @sa mesh::Entity
   * @{
   */
  geometry::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }
  bool operator==(const mesh::Entity& rhs) const override { return this == &rhs; }
  /** @} */
  
  virtual ~Quadrilateral() override = default;
 private:
  size_type index_ = -1;                         // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_; // shape information
  std::array<const Node *,4> corners_;           // nodes = corners of quad
  std::array<const Edge *,4> edges_;             // edges of quad
};    
/** @} */
  
/**
 * @brief classes for topological entities in a 2D hybrid mesh
 * @tparam CODIM the co-dimension of the entity object \f$\in\{0,1,2\}\f$
 *
 * This class template can be used to instantiate the four different topological
 * entities occurring in 2D hybrid meshes: Points, Edges, Triangles, and
 * Quadrilaterals.
 * @note Every `Entity` object owns a smart pointer to an associated geometry
 * object.
 *
 */
template <char CODIM>
// NOLINTNEXTLINE(hicpp-member-init)
class Entity : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Entity() = default;

  Entity(const Entity&) = delete;
  Entity(Entity&&) noexcept = default;
  Entity& operator=(const Entity&) = delete;
  Entity& operator=(Entity&&) noexcept = default;

  // constructor, is called from Mesh
  explicit Entity(
      Mesh* mesh, size_type index,

      /**
       * @brief constructor, is called from MeshFactory
       * @param mesh pointer to global hybrid mesh object
       * @param index index of the entity to be created; will usually be
       * retrieved via the `Index()` method of `Mesh`
       * @param geometry pointer to a geometry object providing the shape of the
       * entity
       * @param sub_entities indices of the sub-entities in the entity arrays of
       * the global mesh
       *
       * @note Note that you need to create a suitable geometry object for the
       * entity before you can initialize the entity object itseld.
       */
      std::unique_ptr<geometry::Geometry>&& geometry,
      std::array<std::vector<size_type>, 2 - CODIM> sub_entities)
      : mesh_(mesh),
        index_(index),
        geometry_(std::move(geometry)),
        sub_entities_(std::move(sub_entities)) {}

  char Codim() const override { return CODIM; }

  base::RandomAccessRange<const mesh::Entity> SubEntities(
      char rel_codim) const override;

  geometry::Geometry* Geometry() const override { return geometry_.get(); }

  base::RefEl RefEl() const override {
    switch (CODIM) {
      case 0:
        return sub_entities_[0].size() == 3 ? base::RefEl::kTria()
                                            : base::RefEl::kQuad();
      case 1:
        return base::RefEl::kSegment();
      case 2:
        return base::RefEl::kPoint();
      default:
        LF_VERIFY_MSG(false, "codim out of range.");
    }
  }

  bool operator==(const mesh::Entity& rhs) const override {
    return this == &rhs;
  }

  ~Entity() override = default;

 private:
  Mesh* mesh_ = nullptr;  // pointer to global hybrid 2D mesh object
  size_type index_ = -1;  // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_;
  std::array<std::vector<size_type>, 2 - CODIM> sub_entities_;

  friend class Mesh;
};

}  // namespace lf::mesh::hybrid2d

#include "mesh.h"

namespace lf::mesh::hybrid2d {

template <char CODIM>
base::RandomAccessRange<const mesh::Entity> Entity<CODIM>::SubEntities(
    char rel_codim) const {
  switch (2 - CODIM - rel_codim) {
    case 2:
      // This case is relevant only for CODIM = 0 and codim =0,
      // that is for cells; return ourselves as the only element of the range
      return {this, this + 1};
    case 1:
      // This case is visited, if
      // (i) either the entity is an edge (CODIM = 1)
      return {base::make_DereferenceLambdaRandomAccessIterator(
                  sub_entities_[rel_codim - 1].begin(),
                  [&](auto i) -> const mesh::Entity& {
                    return mesh_->entities1_[*i];
                  }),
              base::make_DereferenceLambdaRandomAccessIterator(
                  sub_entities_[rel_codim - 1].begin(),
                  [&](auto i) -> const mesh::Entity& {
                    return mesh_->entities1_[*i];
                  })};
    case 0:
      return {base::make_DereferenceLambdaRandomAccessIterator(
                  sub_entities_[rel_codim - 1].begin(),
                  [&](auto i) -> const mesh::Entity& {
                    return mesh_->entities2_[*i];
                  }),
              base::make_DereferenceLambdaRandomAccessIterator(
                  sub_entities_[rel_codim - 1].begin(),
                  [&](auto i) -> const mesh::Entity& {
                    return mesh_->entities2_[*i];
                  })};
    default:
      LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}

}  // namespace lf::mesh::hybrid2d

#endif  // __7a3f1903d42141a3b1135e8e5ad72c1c
