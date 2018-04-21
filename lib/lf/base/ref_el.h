/// @file ref_el.h

#ifndef __96e6ff0ee0034f4584fcdfc7e9c53f82
#define __96e6ff0ee0034f4584fcdfc7e9c53f82
#include <exception>
#include <vector>
#include <array>

#include <boost/range.hpp>
#include <Eigen/Eigen>
#include "lf_assert.h"


namespace lf::base {


/*! 
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
enum class RefElType {
  kPoint,
  //!< @copydoc RefEl::kPoint
  kSegment,
  //!< @copydoc RefEl::kSegment
  kTria,
  //!< @copydoc RefEl::kTria
  kQuad,
  //!< @copydoc RefEl::kQuad
};

/**
 * @brief This is a test class.
 */
class RefEl {
private:
  // Member variable
  RefElType type_;

  // Node coordinates as dynamic eigen matrices:
  static const std::vector<Eigen::VectorXd> ncoords_point_dynamic_;
  static const std::vector<Eigen::VectorXd> ncoords_segment_dynamic_;
  static const std::vector<Eigen::VectorXd> ncoords_tria_dynamic_;
  static const std::vector<Eigen::VectorXd> ncoords_quad_dynamic_;

  // Node coordinates as fixed eigen matrices:
  static const std::vector<Eigen::Matrix<double, 0, 1>> ncoords_point_static_;
  static const std::vector<Eigen::Matrix<double, 1, 1>> ncoords_segment_static_;
  static const std::vector<Eigen::Vector2d> ncoords_tria_static_;
  static const std::vector<Eigen::Vector2d> ncoords_quad_static_;

  // subSubEntities, used by SubSubEntity2SubEntity
  static constexpr std::array<std::array<char,2>,3> sub_sub_entity_index_tria_ = {{{0,1}, {1,2}, {2,3}}};
  static constexpr std::array<std::array<char,2>,4> sub_sub_entity_index_quad_ = {{{0,1},{1,2},{2,3},{3,0}}};

  // Some utility methods that are needed to deduce compile time return types:
  static constexpr char DimensionImpl(RefElType type) {
    switch (type) {
      case RefElType::kPoint:
        return 0;
      case RefElType::kSegment:
        return 1;
      case RefElType::kTria:
        return 2;
      case RefElType::kQuad:
        return 2;
      default:
        throw std::
          exception("RefEl::Dimension() not implemented for this RefEl type.");
    }
  }

public:

  /**
   * @brief Type of the node coordinate iterator that is returned from
   * the static version of NodeCoords()
   */
  template <RefElType type>
  using NodeCoordVector = Eigen::Matrix<double, DimensionImpl(type), 1>;


  /**
   * @brief Create a RefEl from a lf::base::RefElType enum.
   * @param type The type of the Reference Element
   */
  constexpr RefEl(RefElType type) : type_(type) {
  }

  constexpr RefEl(const RefEl&) = default;

  /**
   * @brief The (0-dimensional) reference point.
   */
  static const RefEl kPoint;

  /**
   * @brief The (1-dimensional) reference segment.
   *
   * ### Node numbering with (1D) node coordinates
   * @image html segment.png
   */
  static const RefEl kSegment;


  /**
   * @brief The reference triangle
   *
   * ### Node numbering with (2D) node coordinates and segment orientation.
   * @image html tria.png
   */
  static const RefEl kTria;


  /**
   * @brief The reference quadrilateral
   *
   * ### Node numbering with (2D) node coordinates and segment orientation
   * @image html quad.png
   */
  static const RefEl kQuad;


  /**
   * @brief Return the dimension of this reference element.
   *
   * - 0 for a RefEl::kPoint
   * - 1 for a RefEl::kSegment
   * - 2 for a RefEl::kTria
   * - 2 for a RefEl::kQuad
   */
  constexpr char Dimension() const {
    return DimensionImpl(type_);
  }


  /**
   * @brief The number of nodes of this reference element.
   */
  constexpr char NumNodes() const {
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
        throw std::
          exception("RefEl::NumNodes() not implemented for this RefEl type.");
    }
  }


  /**
   * @brief Get the coordinates of the nodes of this reference element.
   * @return Returns a `std::vector` with `NumNodes()` elements.
   * Every element is a `Eigen::VectorXd` with `Dimension()` rows.
   *   
   * @remark This method is not optimal from a performance point of view
   * because the vector and the Eigen matrices are allocated on the heap.
   * If the type of the RefEl is known at compile time, use that
   * static NodeCoords() function instead.
   * 
   * @snippet ref_el.cc nodeCoordStatic
   */
  const std::vector<Eigen::VectorXd>& NodeCoords() const {
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
      LF_VERIFY_MSG(false,
          "RefEl::NodeCoords() not implemented for this RefEl type.");
    }
  }

  /**
   * @brief Get the coordinates of the nodes of a reference element.
   * @tparam type The `RefElType` of the reference element for which the node
   * coordinates should be returned, see usage below.
   * @return Returns a `std::vector` with `NumNodes()` elements. Every
   * Element is a fixed-size column vector (e.g. `Eigen::Matrix<double,1,1>` 
   * for a segment or `Eigen::Matrix<double,2,1>` for a triangle/quadrilateral)
   * 
   * #### Usage example:
   * @snippet ref_el.cc nodeCoordStatic
   * 
   * @remark This function template can only be called if the type of the
   * reference element is known at compile time. The advantage of this method
   * over NodeCoords() const is, that it returns a vector of fixed-size 
   * vectors that are allocated on the stack and are much easier copied around.
   */
  template <RefElType type>
  static const std::vector<NodeCoordVector<type>>& NodeCoords() {
    if constexpr (type == RefElType::kPoint) {
      return ncoords_point_static_;
    } else if constexpr (type == RefElType::kSegment) {
      return ncoords_segment_static_;
    } else if constexpr (type == RefElType::kTria) {
      return ncoords_tria_static_;
    } else if constexpr (type == RefElType::kQuad) {
      return ncoords_quad_static_;
    } else {
      LF_VERIFY_MSG(false,
        "RefEl::NodeCoords<>() not implemented for this RefEl type.");
    }
  }

  /**
   * @brief Get the number of sub-entities of the given codimension.
   * @param sub_codim The codimension of the subEntities that should be counted.
   * 
   * #### Examples
   * - a segment has two points as `codim=1` sub-entities, therefore 
   *   `RefEl::kSegment.NumSubEntities(1) == 2`
   * - a Triangle has three subentities of `codim=1` (all Segments), therefore
   *  `RefEl::kTria.NumSubEntities(1) == 3`
   * - a Triangle has three subEntities of `codim=2` (all Points), therefore
   *   `RefEl::kTria.NumSubEntities(2) == 3`
   */
  constexpr char NumSubEntities(char sub_codim) const {
    LF_ASSERT_MSG_CONSTEXPR(sub_codim >= 0, "sub_codim is negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_codim <= Dimension(), "sub_codim > Dimension()"
    );
    if (sub_codim == 0) return 1;
    switch (type_) {
      case RefElType::kSegment:
        return 2; // sub_codim=1
      case RefElType::kTria:
        return 3; // sub_codim=1,2
      case RefElType::kQuad:
        return 4; // sub_codim=1,2
      default:
      LF_VERIFY_MSG(false,
          "RefEl::NumSubEntities() not implemented for this RefElType.");
    }
  }

  constexpr RefEl SubType(char sub_codim, char sub_index) const {
    LF_ASSERT_MSG_CONSTEXPR(sub_codim >= 0, "sub_codim is negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_codim <= Dimension(), "sub_codim > Dimension()"
    );
    LF_ASSERT_MSG_CONSTEXPR(sub_index >= 0, "sub_index is negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_index < NumSubEntities(sub_codim),
      "sub_index >= NumSubEntities");

    if (sub_codim == 0) return *this;
    if (sub_codim == Dimension()) return RefElType::kPoint;
    if (Dimension() - sub_codim == 1) return RefElType::kSegment;
    LF_VERIFY_MSG(false, "This code should never be reached.");
  }

  constexpr char SubSubEntity2SubEntity(char sub_codim, char sub_index,
                                        char sub_sub_codim,
                                        char sub_sub_index) const {
    LF_ASSERT_MSG_CONSTEXPR(sub_codim >= 0, "sub_codim negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_codim <= Dimension(), "sub_codim > Dimension");
    LF_ASSERT_MSG_CONSTEXPR(sub_index >= 0, "sub_index negative");
    LF_ASSERT_MSG_CONSTEXPR(sub_index <= NumSubEntities(sub_codim),
      "sub_index >= NumSubEntities");
    LF_ASSERT_MSG_CONSTEXPR(sub_sub_codim >= 0, "sub_sub_codim negative.");
    LF_ASSERT_MSG_CONSTEXPR(sub_sub_codim <= Dimension() - sub_codim,
      "subSubCodim out of bounds.");
    LF_ASSERT_MSG_CONSTEXPR(sub_sub_index >= 0, "sub_sub_index negative.");
    LF_ASSERT_MSG_CONSTEXPR(sub_sub_index < SubType(sub_codim, sub_index).
      NumSubEntities(sub_sub_codim), "sub_sub_index out of bounds.");

    if(type_ == RefElType::kPoint) return 0;
    if(sub_codim == 0) return sub_sub_index;
    if(sub_codim == Dimension()) return sub_index;

    // from here on, it must be a segment
    switch (type_)
    {
    case RefElType::kTria:
      return sub_sub_entity_index_tria_[sub_index][sub_sub_index];
    case RefElType::kQuad:
      return sub_sub_entity_index_quad_[sub_index][sub_sub_index];
    }

    LF_VERIFY_MSG(false, "This code should never be reached.");
  }

  /**
   * @brief Conversion operator, converts this RefEl to a lf::base::RefElType
   * enum.
   *
   * #### Usage example
   * @snippet ref_el.cc convert_to_enum
   */
  constexpr operator RefElType() const {
    return type_;
  }


};

constexpr RefEl RefEl::kPoint = RefEl(RefElType::kPoint);
constexpr RefEl RefEl::kSegment = RefEl(RefElType::kSegment);
constexpr RefEl RefEl::kTria = RefEl(RefElType::kTria);
constexpr RefEl RefEl::kQuad = RefEl(RefElType::kQuad);

}

#endif // __96e6ff0ee0034f4584fcdfc7e9c53f82
