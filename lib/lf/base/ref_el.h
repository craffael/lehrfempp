/// @file ref_el.h

#ifndef __96e6ff0ee0034f4584fcdfc7e9c53f82
#define __96e6ff0ee0034f4584fcdfc7e9c53f82

#include <array>
#include <vector>

//#include <boost/range.hpp>

#include "base.h"
#include "lf_assert.h"

#include <Eigen/Eigen>

namespace lf::base {

/**
 * @brief An enum that defines all possible RefEl types.
 *
 * This enum is only rarely used direcly because there is a one-to-one relation
 * between every enum value and an instance of the lf::base::RefEl class:
 * @snippet ref_el.cc oneToOneRelation
 *
 * Also the enum representation is convertible into a lf::base::RefEl instance
 * and back:
 * @snippet ref_el.cc enumConversion
 *
 */
enum class RefElType : unsigned char {
  kPoint = 1,
  //!< @copydoc RefEl::kPoint()
  kSegment = 2,
  //!< @copydoc RefEl::kSegment()
  kTria = 3,
  //!< @copydoc RefEl::kTria()
  kQuad = 4,
  //!< @copydoc RefEl::kQuad()
};

namespace internal {
// Some utility methods that are needed to deduce compile time return types:
constexpr dim_t DimensionImpl(RefElType type) {
  switch (type) {
    case RefElType::kPoint:
      return 0;
    case RefElType::kSegment:
      return 1;
    case RefElType::kTria:
    case RefElType::kQuad:
      return 2;
    default:
      throw std::runtime_error(
          "RefEl::Dimension() not implemented for this RefEl type.");
  }
}
}  // namespace internal

/** @class RefEl lf/base/base.h
 * @brief Represents a reference element with all its properties.
 *
 * Every entity of a mesh in LehrFEM++ is the image of a reference element
 * under an entity-specific (smooth) transformation (which is described by the
 * lf::mesh::Geometry class). This transformation describes the shape of the
 * actual entity, but also the algebraic relations between its sub-entities.
 * For details consult [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{sec:trftech. @lref{def:parfe}
 *
 * On a more fundamental level objects of this class distinguish the
 * **topological type** of a mesh entity, see [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{par:loctop}.
 *
 * Summing up, the reference element serves two main purposes:
 * - it defines the _local topology_ of the mesh in terms of incidence relations
 *   of the entities. For instance the reference element defines the number of
 *   edges of a cell.
 * - it fixes the parameter domain for the mapping describing the _shape_ of the
 *   entitiy.
 *
 * There is a fixed number of reference elements. This class has a static
 * member field for every type of reference element:
 * - `RefEl::kPoint()` is the Reference element of every point/node in a mesh.
 *   The point itself doesn't have any sub-entities.
 * - `RefEl::kSegment()` is the reference element of every edge in a mesh.
 *   It connects two points with each other.
 * - `RefEl::kTria()` is the reference element of every triangular element in
 * the mesh. It has three segments (codim=1) and three points (codim=2) as
 *   sub-entities.
 * - `RefEl::kQuad()` is the reference element of every quadrilateral element in
 *   the mesh. It has four segments (codim=1) and four points (codim=2) as
 *   sub-entities.
 *
 *
 * #### Usage of this class
 * - You can create arbitrary many instances of this class, but every instance
 *   is equal to one of the four types of reference elements that are exposed
 *   as static member fields (see above).
 * - Instances of this class are copyable, moveable, assignable and have a
 *   number of member functions that give information about the reference
 *   element they represent.
 * - This class is very lightweight, in fact `sizeof(RefEl) == sizeof(char)`.
 *   It can be copied around as needed.
 *
 * @snippet ref_el.cc refElUsage
 */
class RefEl {
  // Node coordinates as dynamic eigen matrices:
  static const Eigen::MatrixXd ncoords_point_dynamic_;
  static const Eigen::MatrixXd ncoords_segment_dynamic_;
  static const Eigen::MatrixXd ncoords_tria_dynamic_;
  static const Eigen::MatrixXd ncoords_quad_dynamic_;

  // Node coordinates as fixed eigen matrices:
  static const std::vector<Eigen::Matrix<double, 0, 1>> ncoords_point_static_;
  static const std::vector<Eigen::Matrix<double, 1, 1>> ncoords_segment_static_;
  static const std::vector<Eigen::Vector2d> ncoords_tria_static_;
  static const std::vector<Eigen::Vector2d> ncoords_quad_static_;

  // Member variable
  RefElType type_;

  // subSubEntities, used by SubSubEntity2SubEntity
  static constexpr std::array<std::array<base::dim_t, 2>, 3>
      sub_sub_entity_index_tria_ = {{{0, 1}, {1, 2}, {2, 0}}};
  static constexpr std::array<std::array<base::dim_t, 2>, 4>
      sub_sub_entity_index_quad_ = {{{0, 1}, {1, 2}, {2, 3}, {3, 0}}};

 public:
  using dim_t = unsigned int;
  /**
   * @brief Type of the node coordinate iterator that is returned from
   * NodeCoords()
   */
  template <RefElType type>
  using NodeCoordVector =
      Eigen::Matrix<double, internal::DimensionImpl(type), 1>;

  /**
   * @brief Returns the (0-dimensional) reference point.
   */
  static constexpr RefEl kPoint() { return RefEl(RefElType::kPoint); }

  /**
   * @brief Returns the (1-dimensional) reference segment.
   *
   * #### Node numbering with (1D) node coordinates
   *
   * @image html segment_num_edg.svg
   */
  static constexpr RefEl kSegment() { return RefEl(RefElType::kSegment); }

  /**
   * @brief Returns the reference triangle
   *
   * #### Node numbering with (2D) node coordinates and segment orientation.
   * @image html tria_num_edg.svg
   */
  static constexpr RefEl kTria() { return RefEl(RefElType::kTria); }

  /**
   * @brief Returns the reference quadrilateral
   *
   * #### Node numbering with (2D) node coordinates and segment orientation
   * @image html quad_num_edg.svg
   */
  static constexpr RefEl kQuad() { return RefEl(RefElType::kQuad); }

  /**
   * @brief Create a RefEl from a lf::base::RefElType enum.
   * @param type The type of the Reference Element
   */
  explicit constexpr RefEl(RefElType type) noexcept : type_(type) {}

  /** @brief Default copy constructor */
  constexpr RefEl(const RefEl&) = default;

  /** @brief Default move constructor */
  constexpr RefEl(RefEl&&) = default;

  /** @brief Default copy assignment operator */
  // NOLINTNEXTLINE(modernize-use-equals-default,hicpp-use-equals-default,cert-oop54-cpp)
  constexpr RefEl& operator=(const RefEl& rhs) {
    type_ = rhs.type_;
    return *this;
  }

  /** @brief Default move assignment operator */
  constexpr RefEl& operator=(RefEl&& rhs) noexcept {
    type_ = rhs.type_;
    return *this;
  }

  /**
   * @brief Return the dimension of this reference element.
   *
   * - 0 for a RefEl::kPoint()
   * - 1 for a RefEl::kSegment()
   * - 2 for a RefEl::kTria()
   * - 2 for a RefEl::kQuad()
   */
  [[nodiscard]] constexpr dim_t Dimension() const {
    return internal::DimensionImpl(type_);
  }

  /**
   * @brief The number of nodes of this reference element.
   *
   * @remark This is a shortcut for calling `NumSubEntities(Dimension())`
   */
  [[nodiscard]] constexpr size_type NumNodes() const {
    switch (type_) {
      case RefElType::kPoint:
        return 1;
      case RefElType::kSegment:
        return 2;
      case RefElType::kTria:
        return 3;
      case RefElType::kQuad:
        return 4;
      default:
        throw std::runtime_error(
            "RefEl::NumNodes() not implemented for this RefEl type.");
    }
  }

  /**
   * @brief Get the coordinates of the nodes of this reference element.
   * @return a matrix of size `Dimension() x NumNodes()` that contains the
   *         the coordinates of every node as a column vector.
   *
   * @remark This method is not optimal from a performance point of view
   * because the matrix is allocated on the heap.
   * If the type of the RefEl is known at compile time, use the
   * static NodeCoords() function instead.
   *
   * @snippet ref_el.cc nodeCoordStatic
   */
  [[nodiscard]] const Eigen::MatrixXd& NodeCoords() const {
    switch (type_) {
      case RefElType::kPoint:
        return ncoords_point_dynamic_;
      case RefElType::kSegment:
        return ncoords_segment_dynamic_;
      case RefElType::kTria:
        return ncoords_tria_dynamic_;
      case RefElType::kQuad:
        return ncoords_quad_dynamic_;
      default:
        LF_VERIFY_MSG(
            false, "RefEl::NodeCoords() not implemented for this RefEl type.");
    }
  }

  /**
   * @brief Get the coordinates of the nodes of a reference element.
   * @tparam type The `RefElType` of the reference element for which the node
   * coordinates should be returned, see usage below.
   * @return a `std::vector` with `NumNodes()` elements. Every
   * Element is a fixed-size column vector (e.g. `Eigen::Matrix<double,1,1>`
   * for a segment or `Eigen::Matrix<double,2,1>` for a triangle/quadrilateral)
   *
   * #### Usage example:
   * @snippet ref_el.cc nodeCoordStatic
   *
   * @remark This function template can only be called if the type of the
   * reference element is known at compile time. The advantage of this method
   * over NodeCoords() const is, that it returns a vector of fixed-size
   * vectors that are allocated on the stack.
   */
  template <RefElType type>
  static const std::vector<NodeCoordVector<type>>& NodeCoords() {
    // NOLINTNEXTLINE
    if constexpr (type == RefElType::kPoint) {
      return ncoords_point_static_;
    }
    // NOLINTNEXTLINE
    if constexpr (type == RefElType::kSegment) {
      return ncoords_segment_static_;
    }
    // NOLINTNEXTLINE
    if constexpr (type == RefElType::kTria) {
      return ncoords_tria_static_;
    }
    // NOLINTNEXTLINE
    if constexpr (type == RefElType::kQuad) {
      return ncoords_quad_static_;
    }
    LF_VERIFY_MSG(false,
                  "RefEl::NodeCoords<>() not implemented for this RefEl type.");
  }

  /**
   * @brief Get the number of sub-entities of this RefEl with the given
   *        codimension.
   * @param sub_codim The codimension of the subEntities that should be counted.
   *
   * #### Examples
   * - a segment has two points as `codim=1` sub-entities, therefore
   *   `RefEl::kSegment().NumSubEntities(1) == 2`
   * - a Triangle has three subentities of `codim=1` (all Segments), therefore
   *  `RefEl::kTria().NumSubEntities(1) == 3`
   * - a Triangle has three subEntities of `codim=2` (all Points), therefore
   *   `RefEl::kTria().NumSubEntities(2) == 3`
   */
  [[nodiscard]] constexpr size_type NumSubEntities(dim_t sub_codim) const {
    LF_ASSERT_MSG_CONSTEXPR(sub_codim >= 0, "sub_codim is negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_codim <= Dimension(),
                            "sub_codim > Dimension()");
    if (sub_codim == 0) {
      return 1;
    }
    switch (type_) {
      case RefElType::kSegment:
        return 2;  // sub_codim=1
      case RefElType::kTria:
        return 3;  // sub_codim=1,2
      case RefElType::kQuad:
        return 4;  // sub_codim=1,2
      default:
        LF_ASSERT_MSG_CONSTEXPR(
            false,
            "RefEl::NumSubEntities() not implemented for this RefElType.");
    }
    return 0;  // prevent warnings from compilers
  }

  /**
   * @brief Return the RefEl of the sub-entity with codim `sub_codim`
   * and index `sub_index`.
   * @param sub_codim The codimension of the sub-entity (w.r.t. `Dimension()`).
   *        Should be `<= Dimension()`.
   * @param sub_index The zero-based index of the sub-entity.
   *        `sub_index` should be smaller than `NumSumEntities(sub_codim)`
   *
   * #### Examples:
   * - A triangle has three codim=2 entities which are all points, therefore
   *   `RefEl::kTria().SubType(2,i) == RefEl::kPoint()` for i=0,1,2.
   * - A quadrilateral has four codim=1 entities which are all segments,
   *   therefore `RefEl::kQuad().SubType(1,i) == RefEl::kSegment()` for
   * i=0,1,2,3.
   * - The codim=0 subEntity of a triangle is the triangle itself, therefore
   *   `RefEl::kTria().SubType(0,0) == RefEl::kTria()`.
   *
   *
   * @see NumSubEntities() const to get the number of sub entities of a given
   *      codimension
   */
  [[nodiscard]] constexpr RefEl SubType(dim_t sub_codim,
                                        dim_t sub_index) const {
    LF_ASSERT_MSG_CONSTEXPR(sub_codim >= 0, "sub_codim is negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_codim <= Dimension(),
                            "sub_codim > Dimension()");
    LF_ASSERT_MSG_CONSTEXPR(sub_index >= 0, "sub_index is negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_index < NumSubEntities(sub_codim),
                            "sub_index >= NumSubEntities");

    if (sub_codim == 0) {
      return *this;
    }
    if (sub_codim == Dimension()) {
      return kPoint();
    }
    if (Dimension() - sub_codim == 1) {
      return kSegment();
    } else {  // NOLINT(readability-else-after-return)
      LF_ASSERT_MSG_CONSTEXPR(false, "This code should never be reached.");
    }

    return kPoint();  // prevent warnings from compiler
  }

  /**
   * @brief Identifies sub-entities of sub-entities (so-called sub-sub-entities)
   *        with sub-entities.
   *
   * @param sub_codim The codimension of the sub-entity.
   * @param sub_index  The zero-based index of the sub-entity.
   * @param sub_rel_codim The codimension of the sub-sub-entity w.r.t. the
   *        sub-entity identified by `sub_codim` and `sub_index`.
   * @param sub_rel_index The index of the sub-sub-entity w.r.t. the
   *        sub-entity identified by `sub_codim` and `sub_index`.
   * @return The index of the sub-sub-entity w.r.t. this RefEl.
   *
   * @note the argument `sub_rel_codim` is a **relative co-dimension**.
   *
   * @note the argument `sub_rel_index` is a local index in the sub-entity.
   *
   * #### Conventions
   * - For a triangle (the arrow indicates the induced orientation of the edge)
   *   - edge 0 connects vertices 0 -> 1
   *   - edge 1 connects vertices 1 -> 2
   *   - edge 2 connects vertices 2 -> 0
   * - For a quadrilateral
   *   - edge 0 connects vertices 0 -> 1
   *   - edge 1 connects vertices 1 -> 2
   *   - edge 2 connects vertices 2 -> 3
   *   - edge 3 connects vertices 3 -> 0
   *
   * #### Examples
   * - The sub-entity of a `RefEl::kTria()` with `sub_codim=1`, `sub_index=1`
   *   is a `RefEl::kSegment()` that connects node 1 with node 2. The
   *   (sub-)sub-entity of this segment with codim `sub_rel_codim=1` and
   *   sub-index `sub_rel_index=0` (both w.r.t. to the segment) is the first
   * point of the segment, i.e. node 1. Therefore
   * `SubSubEntity2SubEntity(1,1,1,0) == 1`
   * - Similarly, for `sub_rel_index=1`:
   *   `SubSubEntity2SubEntity(1,1,1,1) == 2`
   */
  [[nodiscard]] constexpr sub_idx_t SubSubEntity2SubEntity(
      dim_t sub_codim, sub_idx_t sub_index, dim_t sub_rel_codim,
      sub_idx_t sub_rel_index) const {
    LF_ASSERT_MSG_CONSTEXPR(sub_codim >= 0, "sub_codim negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_codim <= Dimension(), "sub_codim > Dimension");
    LF_ASSERT_MSG_CONSTEXPR(sub_index >= 0, "sub_index negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_index <= NumSubEntities(sub_codim),
                            "sub_index >= NumSubEntities");
    LF_ASSERT_MSG_CONSTEXPR(sub_rel_codim >= 0, "sub_rel_codim negative.");
    LF_ASSERT_MSG_CONSTEXPR(sub_rel_codim <= Dimension() - sub_codim,
                            "subSubCodim out of bounds.");
    LF_ASSERT_MSG_CONSTEXPR(sub_rel_index >= 0, "sub_rel_index negative.");
    LF_ASSERT_MSG_CONSTEXPR(
        sub_rel_index <
            SubType(sub_codim, sub_index).NumSubEntities(sub_rel_codim),
        "sub_sub_index out of bounds.");

    if (type_ == RefElType::kPoint) {
      return 0;
    }
    if (sub_codim == 0) {
      return sub_rel_index;
    }
    if (sub_codim == Dimension()) {
      return sub_index;
    }

    // from here on, it must be a segment
    switch (type_) {
      case RefElType::kTria:
        return sub_sub_entity_index_tria_[sub_index][sub_rel_index];
      case RefElType::kQuad:
        return sub_sub_entity_index_quad_[sub_index][sub_rel_index];
      default:
        LF_ASSERT_MSG_CONSTEXPR(false, "This code should never be reached.");
    }

    return 0;  // Prevent warnings from compiler...
  }

  /**
   * @brief Return a string representation of this Reference element
   *
   * This string is supposed to described  the _topological type_ of the entity:
   * NODE, EDGE (2 nodes), TRIA (3-node triangle), QUAD (4-node quadrilateral)
   */
  [[nodiscard]] std::string ToString() const {
    switch (type_) {
      case RefElType::kPoint:
        return "NODE";
      case RefElType::kSegment:
        return "EDGE";
      case RefElType::kTria:
        return "TRIA";
      case RefElType::kQuad:
        return "QUAD";
      default:
        LF_VERIFY_MSG(false, "ToString() not implemented for this RefElType");
    }
  }

  /**
   * @brief Conversion operator, converts this RefEl to a lf::base::RefElType
   * enum.
   *
   * #### Usage example
   * @snippet ref_el.cc convert_to_enum
   */
  // NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions)
  [[nodiscard]] constexpr operator RefElType() const { return type_; }

  /**
   * @brief Return a unique id for this reference element.
   *
   * The id numbers are contiguous over the known entity types starting from 0.
   * Thus, this identification number can be used as an array index, if one
   * wants to store information for (all) type of entities.
   *
   * #### Usage example
   * @snippet ref_el.cc id
   */
  // NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions)
  [[nodiscard]] constexpr unsigned int Id() const {
    return static_cast<unsigned int>(type_);
  }

  ~RefEl() = default;

};  // class RefEl

// Declare print function
/**
 * @brief Diagnostic output operator. Prints info about a reference element.
 * @param ref_el The reference element to print info about
 * @param o The stream to which this function should output
 * @param output_ctrl Determines the level of detail with which the output is
 * printed to `o`
 *
 * #### Output levels
 * - `output_ctrl` = 0: Type of reference element, dimension and number of
 * nodes
 * - `output_ctrl` > 0: The above and number of subentities and their
 * types for each codimension
 * - `output_ctrl` > 10: The above and type of subentity for each
 * subentities in each codimension.
 * - `output_ctrl` > 20: The above and coordinates of the points of the
 * reference element
 *
 * @relates RefEl
 */
void PrintInfo(std::ostream& o, const RefEl& ref_el, int output_ctrl = 0);

/**
 * @brief Operator overload to print a `RefEl` to a stream, such as `std::cout`
 * @param stream The stream to which this function should output
 * @param ref_el The reference element to write to `stream`.
 * @return The stream itself.
 *
 * @note This function directly calls `RefEl::ToString()`.
 *
 * #### Usage example
 * @snippet ref_el.cc streamOutput
 */
inline std::ostream& operator<<(std::ostream& stream, const RefEl& ref_el) {
  return stream << ref_el.ToString();
}

}  // namespace lf::base

#endif  // __96e6ff0ee0034f4584fcdfc7e9c53f82
