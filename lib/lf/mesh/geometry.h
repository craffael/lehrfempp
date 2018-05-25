
#ifndef __7ed6b0d4d9244155819c464fc4eb9bbb
#define __7ed6b0d4d9244155819c464fc4eb9bbb

#include <lf/base/ref_el.h>
#include <Eigen/Eigen>

namespace lf::mesh {

class Geometry {
 public:

  virtual char DimFrom() const = 0;

  virtual char DimTo() const = 0;

  virtual base::RefEl RefEl() const = 0;

  virtual Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const = 0;

  virtual Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const = 0;

  virtual Eigen::VectorXd IntegrationElement(const Eigen::MatrixXd& local) const = 0;
};

}  // namespace lf::mesh

#endif  // __7ed6b0d4d9244155819c464fc4eb9bbb
