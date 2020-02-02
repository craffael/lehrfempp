#ifndef __7ed6b0d4d9244155819c464fc4eb9bbb
#define __7ed6b0d4d9244155819c464fc4eb9bbb

#include <lf/base/base.h>
#include <Eigen/Eigen>
#include <memory>
#include "refinement_pattern.h"

namespace lf::geometry {

/** @brief **Interface class** for shape information on a mesh cell in the
 * spirit of _parametric finite element methods_
 *
 * This abstract base class is supposed to provide complete information about
 * the mapping \f$\Phi:\widehat{K}\to K\f$ that takes the reference cell to a
 * particular cell. The reference cell is defined through the topological type
 * of a cell, see @ref lf::base::RefEl.
 */

class Geometry {
 protected:
  Geometry() = default;
  Geometry(const Geometry&) = default;
  Geometry(Geometry&&) = default;
  Geometry& operator=(const Geometry&) = default;
  Geometry& operator=(Geometry&&) = default;

 public:
  using dim_t = base::RefEl::dim_t;

  /**
   * @brief Dimension of the domain of this mapping.
   */
  [[nodiscard]] virtual dim_t DimLocal() const = 0;

  /**
   * @brief Dimension of the image of this mapping.
   */
  [[nodiscard]] virtual dim_t DimGlobal() const = 0;

  /**
   * @brief The Reference element that defines the domain of this mapping.
   */
  [[nodiscard]] virtual base::RefEl RefEl() const = 0;

  /**
   * @brief Map a number of points in local coordinates into the global
   *        coordinate system.
   * @param local A Matrix of size `DimLocal() x numPoints` that contains
   *              in its columns the coordinates of the points at which the
   *              mapping function should be evaluated.
   * @return A Matrix of size `DimGlobal() x numPoints` that contains the mapped
   *         points as column vectors. Here `numPoints` is the number of columns
   *         of the matrix passed in the `local` argument.
   *
   * This method provides a complete description of the shape of an entity
   * through a parameterization over the corresponding reference element =
   * parameter domain. The method takes as arguments a number of coordinate
   * vectors of points in the reference element. For the sake of efficiency,
   * these coordinate vectors are passed as the columns of a dynamic matrix type
   * as supplied by Eigen.
   *
   * For instance, this method is used in @ref lf::geometry::Corners()
   * ~~~
   inline Eigen::MatrixXd Corners(const Geometry& geo) {
   return geo.Global(geo.RefEl().NodeCoords()); }
   * ~~~
   * Additional explanations in [Lecture
   * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
   * @lref{par:dunetrf}.
   *
   */
  [[nodiscard]] virtual Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const = 0;

  /**
   * @brief Evaluate the jacobian of the mapping simultaneously at `numPoints`
   *        points.
   * @param local A Matrix of size `DimLocal x numPoints` that contains the
   *              evaluation points as column vectors
   * @return A Matrix of size `DimGlobal() x (DimLocal() * numPoints)` that
   * contains the jacobians at the evaluation points.
   *
   * This method allows access to the derivative of the parametrization mapping
   * in a number of points, passed as the columns of a dynamic matrix. The
   * derivative of the parametrization in a point is a Jacobian matrix of size
   * `DimGlobal() x DimLocal()'. For the sake of efficiency, these matrices are
   * stacked horizontally and returned as one big dynamic matrix. Use Eigen's
   * `block()' method of `Eigen::MatrixXd` to extract the individual Jacobians
   * from the returned matrix.
   */
  [[nodiscard]] virtual Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const = 0;

  /**
   * @brief Evaluate the Jacobian * Inverse Gramian (\f$ J (J^T J)^{-1}\f$)
   *        simulatanesouly at `numPoints`.
   * @param local A Matrix of size `DimLocal() x numPoints` that contains the
   *              evaluation points as column vectors.
   * @return A Matrix of size `DimGlobal() x (DimLocal() * numPoints)` that
   *         contains the Jacobian multiplied with the Inverse Gramian
   *         ( \f$ J (J^T J)^{-1}\f$) at every evaluation point.
   *
   * @note If `DimLocal() == DimGlobal()` then \f$ J (J^T J)^{-1} = J^{-T} \f$,
   *       i.e. this method returns the inverse of the transposed jacobian.
   *
   * ### Example for recovering a single transposed inverse Jacobian.
   *
   * If both dimensions agree and have the value D, then the method
   * returns the transposed of the inverse Jacobians of the transformation at
   * the passed points. These are square DxD matrices.
   *
   * To retrieve the j-th inverse transposed Jacobian from the returned matrix,
   * use the `block` methdod of `Eigen`:
   * ~~~
   * JacobianInverseGramian(local).block(0,i*D,D,D)
   * ~~~
   */
  [[nodiscard]] virtual Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const = 0;

  /**
   * @brief The integration element (factor appearing in integral transformation
   *        formula, see below) at number of evaluation points (specified in
   *        local coordinates).
   * @param local A Matrix of size `DimLocal() x numPoints` that contains the
   *         evaluation points (in local = reference coordinates) as column
   *         vectors.
   * @return A Vector of size `numPoints x 1` that contains the integration
   *         elements at every evaluation point.
   *
   * For a transformation \f$ \Phi : K \mapsto R^{\text{DimGlobal}}\f$ with
   * Jacobian \f$ D\Phi : K \mapsto R^{\text{DimGlobal} \times \text{DimLocal}}
   * \f$ the integration element \f$ g \f$ at point \f$ \xi \in K \f$ is defined
   * by
   * \f[
       g(\xi) := \sqrt{\mathrm{det}\left|D\Phi^T(\xi) D\Phi(\xi) \right|}
   * \f]
   *
   * More information also related to the use of lovcal quadrature rules is
   * given in [Lecture
   * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
   * @lref{par:tranquil}.
   */
  [[nodiscard]] virtual Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const = 0;

  /**
   * @brief **Construct** a new Geometry() object that describes the geometry of
   *        the `i`-th sub-entity with codimension=`codim`
   * @param codim The codimension of the sub-entity (w.r.t. `DimLocal()`)
   * @param i The zero-based index of the sub-entity.
   * @return A new Geometry object that describes the geometry of the specified
   *         sub-entity.
   *
   * Let \f$ \mathbf{\Phi} : K \mapsto \mathbb{R}^\text{DimGlobal} \f$ be the
   * mapping of this Geometry object and let \f$ \mathbf{\xi} :
   * \mathbb{R}^{\text{DimLocal}-codim} \mapsto K\f$ be the first-order mapping
   * that maps the reference element `RefEl().SubType(codim,i)` to the `i`-th
   * sub-entity of `RefEl()`. I.e. for every node \f$ \mathbf{n_j} \f$ of
   * `RefEl().SubType(codim,i)` it holds that \f$ \mathbf{\xi}(\mathbf{n_j}) =
   * \f$ `RefEl().NodeCoords(RefEl().SubSubEntity2SubEntity(codim, i,
   * DimLocal()-codim, j))`.
   *
   * Then the geometry element returned by this method describes exactly the
   * mapping \f$ \mathbf{\Phi} \circ \mathbf{\xi} \f$
   */
  [[nodiscard]] virtual std::unique_ptr<Geometry> SubGeometry(
      dim_t codim, dim_t i) const = 0;

  /**
   * @brief **Generate** geometry objects for child entities created in
   *        the course of refinement.
   *
   * @param ref_pat An object encoding the details of refinement
   * @param codim _relative_ codimension of child entities whose
   *         shape is requested
   * @return an array of unique pointers pointers to geometry objects for child
   *         entities of the specified co-dimension and the given refinement
   *         pattern. The numbering of child entities is a convention.
   *
   * Implementation of local mesh refinement entails splitting of elements into
   * parts ("children"), whose shape will depend on the refinement pattern.
   * This method creates the geometry objects describing the shape of children.
   * The details of subdivisions corresponding to particular refinement patterns
   * are fixed by the method @ref Hybrid2D::RefinementPattern::ChildPolygons()
   * and should be documented there.
   *
   * @sa class lf::geometry::RefinementPattern
   *
   */
  [[nodiscard]] virtual std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat, lf::base::dim_t codim) const = 0;

  /**
   * @brief element shape by affine mapping from reference element
   *
   * @return true, if the element is the affine image of a reference element
   *
   * An affine map is a linear mapping plus a translation. For a 2D mesh
   * an entity has an affine geometry, if
   * - a segement is a straight line
   * - a triangle is flat
   * - a quadrilateral is a flat parallelogram.
   *
   * @note The Jacobian/Gramian for an affine element are _constant_.
   */
  [[nodiscard]] virtual bool isAffine() const { return true; }

  /**
   * @brief Virtual destructor
   */
  virtual ~Geometry() = default;

  // Output control variable
  /** @brief Output control variable */
  static unsigned int output_ctrl_;

};  // class Geometry

/**
 * @brief Compute the (approximate) volume (area) of a shape
 * @param geometry object
 * @return approximate volume
 *
 * @note the volume can be computed exactly only for planar affine/bilinear (2D)
 * shapes Otherwise this functions returns a one-point quadrature approximation
 */
double Volume(const Geometry& geo);

/**
 * @brief The corners of a shape with piecewise smooth boundary
 * @param geometry object
 * @return the coordinates vectors for the corners of the shape represented by
 * the geometry object packed into the _columns_ of a dxn-matrix, where d is the
 * dimension of physical space, and n stands for the number of corners.
 *
 * @note Only for 2D shapes with all straight edges the positions of the corners
 * will completely define the shape. This is the case with shapes represented by
 * the types
 * @ref lf::geometry::TriaO1 and @ref lf::geometry::QuadO1
 *
 * #### Demonstration code [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{cpp:lfgeo}
 * @snippet geometryuse.cc corners
 *
 * Additional explanations can be found in [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{par:coords}.
 *
 */
inline Eigen::MatrixXd Corners(const Geometry& geo) {
  return geo.Global(geo.RefEl().NodeCoords());
}

}  // namespace lf::geometry

#endif  // __7ed6b0d4d9244155819c464fc4eb9bbb
