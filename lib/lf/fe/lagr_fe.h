#ifndef LF_FE
#define LF_FE
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structures representing simple finite elements
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <lf/assemble/dofhandler.h>

namespace lf::fe {
/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::assemble::ldof_idx_t;
/** Type for vector length/matrix sizes */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;

template <typename SCALAR>
class H1ReferenceFiniteElement {
 protected:
  H1ReferenceFiniteElement() = default;
  H1ReferenceFiniteElement(const H1ReferenceFiniteElement&) = default;
  H1ReferenceFiniteElement(H1ReferenceFiniteElement&&) = default;
  H1ReferenceFiniteElement& operator=(const H1ReferenceFiniteElement&) =
      default;
  H1ReferenceFiniteElement& operator=(H1ReferenceFiniteElement&&) = default;

 public:
  explicit H1ReferenceFiniteElement(lf::base::RefEl ref_el);

  lf::base::RefEl RefEl() const { return ref_el_; }
  virtual dim_t Dimension() const = 0;
  virtual size_type NumRefShapeFunctions() const = 0;
  virtual size_type NumRefShapeFunctions(dim_t codim) const = 0;
  virtual size_type NumRefShapeFunctions(dim_t codim,
                                         sub_idx_t subidx) const = 0;
  virtual std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const = 0;
  virtual std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  GradientsReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const = 0;
  virtual Eigen::MatrixXd InterpolationNodes() const = 0;

  virtual ~H1ReferenceFiniteElement() = default;

 private:
  lf::base::RefEl ref_el_;
};
}  // namespace lf::fe

#endif
