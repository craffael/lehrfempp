/**
 * @file
 * @brief Defines C++20 Concepts for the assembly module
 * @author Raffael Casagrande
 * @date January 2024
 * @copyright MIT License
 */

#ifndef INCG_69d30d4083db4e6382d513cd030df97b
#define INCG_69d30d4083db4e6382d513cd030df97b

#include <lf/mesh/mesh.h>

namespace lf::assemble {

/**
 * @ingroup concepts
 * @headerfile lf/assemble/assemble.h
 * @brief Provides the local element matrix for a mesh entity.
 *
 * @tparam EMP Entity-Matrix-Provider (EMP) Type to be tested for satisfaction
 *         of the concept
 *
 * ## Description
 *
 * An EntityMatrixProvider is an object that provides a local element matrix for
 * a given mesh entity. EntityMatrixProvider's are mostly used together with
 * `lf::assemble::AssembleMatrixLocally` to assemble a sparse system matrix.
 *
 * ## Requirements
 *
 * The type `EMP` satisfies the Concept `EntityMatrixProvider` if given
 * - `emp` an object of type `EMP`
 * - `e` a `const lf::mesh::Entity` object
 * - `SCALAR` is an Eigen supported scalar datatype (e.g. `double`)
 * - `ROWS` is the number of rows of the local element matrix or
 * `Eigen::Dynamic`
 * - `COLS` is the number of columns of the local element matrix or
 * `Eigen::Dynamic`
 *
 * the following expressions are valid:
 *
 * <table>
 * <tr>
 *   <th>expression
 *   <th>return type
 *   <th>semantics
 * <tr>
 *   <td> `emp.isActive(e)`
 *   <td> `bool`
 *   <td> Defines whether the entity `e` is taken into account by the
 *        assembly routine.
 * <tr>
 *   <td> `emp.Eval(e)`
 *   <td> `Eigen::Matrix<SCALAR, ROWS, COLS>`
 *   <td> Returns the local element matrix for mesh entity `e`.
 *        Is only called if `emp.isActive(e)==true`.
 * </table>
 *
 * ### Typical class definition
 *
 * @snippet assembler.cc lflinfeelmat
 *
 * @note returning a matrix object through `Eval()` may sacrifice efficiency
 * because it may entail additional allocation of dynamic memory. However, this
 * design was chosen for the sake of readability of the code and in order to
 * avoid nasty memory access errors that often occur when passing a matrix by
 * reference.
 *
 * ## Usage scenarios
 * The concept of a EntityMatrixProvider is widely used in the `lf::assemble`
 * and `lf::uscalfe` modules:
 * - The function `lf::assemble::AssembleMatrixLocally()` accepts an
 *   `EntityMatrixProvider` which in turn defines how the (sparse) system matrix
 *   is assembled.
 *
 * ## Archetype
 * - \ref lf::assemble::EntityMatrixProviderAT<SCALAR> "EntityMatrixProviderAT".
 *
 * ## See also
 * - The concept \ref entity_vector_provider
 * - \ref assemble_matrix_locally "AssembleMatrixLocally()"
 * - [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{ex:lfemp}
 * - [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{ex:lfbdse1}
 *
 * ## Classes modelling this concept
 * - @ref lf::fe::DiffusionElementMatrixProvider
 * - @ref lf::fe::MassElementMatrixProvider
 * - @ref lf::fe::MassEdgeMatrixProvider
 * - @ref lf::uscalfe::LinearFELaplaceElementMatrix
 * - @ref lf::uscalfe::ReactionDiffusionElementMatrixProvider
 * - @ref lf::uscalfe::MassEdgeMatrixProvider
 */
template <class EMP>
concept EntityMatrixProvider = requires(EMP& emp, const mesh::Entity& e) {
  { emp.isActive(e) } -> std::same_as<bool>;
  { emp.Eval(e) } -> base::EigenMatrix<void, -1, -1>;
};

/**
 * @brief \ref archetypes "Archetype" for the EntityMatrixProvider concept.
 * @tparam SCALAR The scalar type of the matrix entries.
 *
 * @sa archetypes
 */
template <class SCALAR>
struct EntityMatrixProviderAT {
  /// Copy constructor deleted because not part of the concept
  EntityMatrixProviderAT(const EntityMatrixProviderAT&) noexcept = delete;

  /// Move constructor deleted because not part of the concept
  EntityMatrixProviderAT(EntityMatrixProviderAT&&) noexcept = delete;

  /// Copy assignment operator deleted because not part of the concept
  auto operator=(const EntityMatrixProviderAT&) noexcept
      -> EntityMatrixProviderAT& = delete;

  /// Move assignment operator deleted because not part of the concept
  auto operator=(EntityMatrixProviderAT&&) noexcept
      -> EntityMatrixProviderAT& = delete;

  /// Destructor
  ~EntityMatrixProviderAT() noexcept = delete;

  bool isActive(const mesh::Entity& e) { return true; }

  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> Eval(
      const mesh::Entity& e) {
    return {};
  }
};

/**
 * @headerfile lf/assemble/assemble.h
 * @brief Provides the local element vector for a mesh entity.
 *
 * ## Description
 *
 * An EntityVectorProvider is an object that provides a local element vector for
 * a given mesh entity. EntityVectorProvider's are mostly used together with
 * `lf::assemble::AssembleVectorLocally` to assemble a vector.
 *
 * ## Requirements
 *
 * The type `EVP` satisfies the Concept `EntityVectorProvider` if given
 * - `evp` an object of type `EVP`
 * - `e` is a const lf::mesh::Entity object
 * - `SCALAR` is an Eigen supported scalar datatype (e.g. `double`)
 * - `ROWS` is the number of rows of the local element vector or
 *   `Eigen::Dynamic`
 *
 * the following expressions are valid:
 *
 * <table>
 * <tr>
 *   <th>expression
 *   <th>return type
 *   <th>semantics
 * <tr>
 *   <td> `evp.isActive(e)`
 *   <td> `bool`
 *   <td> Defines whether the entity `e` is taken into account by the assembly
 *        routine.
 * <tr>
 *   <td> `evp.Eval(e)`
 *   <td> `Eigen::Matrix<SCALAR, ROWS, 1>`
 *   <td> Returns the local element vector for mesh entity `e`. Is only called
 *        if `evp.isActive(e)==true`.
 * </table>
 *
 * ## Usage scenarios
 * The concept of a EntityVectorProvider is widely used in the `lf::assemble`
 * and `lf::uscalfe` modules:
 * - The function `lf::assemble::AssembleVectorLocally()` accepts an
 * `EntityVectorProvider` which in turn defines the global vector is assembled.
 *
 * ## Archetype
 * - \ref lf::assemble::EntityVectorProviderAT "EntityVectorProviderAT"
 *
 * ## See also
 * - The concept lf::assemble::EntityMatrixProvider
 * - \ref assemble_vector_locally "AssembleVectorLocally()"
 *
 * ## Classes modelling this concept
 * - lf::fe::ScalarLoadElementVectorProvider
 * - lf::fe::ScalarLoadEdgeVectorProvider
 * - lf::uscalfe::LinearFELocalLoadVector
 * - lf::uscalfe::ScalarLoadElementVectorProvider
 * - lf::uscalfe::ScalarLoadEdgeVectorProvider
 * 
 */
template <class EVP>
concept EntityVectorProvider = requires(EVP& evp, const mesh::Entity& e) {
  { evp.isActive(e) } -> std::same_as<bool>;
  { evp.Eval(e) } -> base::EigenMatrix<void, -1, 1>;
};

/**
 * @brief \ref archetypes "Archetype" for the EntityVectorProvider concept.
 * @tparam SCALAR The scalar type of the vector entries.
 *
 * @sa archetypes
 */
template <class SCALAR>
struct EntityVectorProviderAT {
  bool isActive(const mesh::Entity& e);

  Eigen::Vector<SCALAR, Eigen::Dynamic> Eval(const mesh::Entity& e);
};

}  // namespace lf::assemble

#endif  // INCG_69d30d4083db4e6382d513cd030df97b