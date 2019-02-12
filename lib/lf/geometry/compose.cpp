/**
 * @file
 * @brief Implementation of the compose function.
 * @author Raffael Casagrande
 * @date   2018-06-17 06:25:57
 * @copyright MIT License
 */

#include "compose.h"

namespace lf::geometry {

namespace /*anonymous */ {

class CompositeGeometry : public Geometry {
 public:
  using dim_t = base::RefEl::dim_t;
  dim_t DimLocal() const override { return dim_local_; }
  dim_t DimGlobal() const override { return dim_global_; }
  base::RefEl RefEl() const override { return b_->RefEl(); }
  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override {
    return a_->Global(b_->Global(local));
  }
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override {
    auto Da = a_->Jacobian(local);
    auto Db = b_->Jacobian(local);
    Eigen::MatrixXd result(DimGlobal(), local.cols() * dim_local_);
    for (int i = 0; i < local.cols(); ++i) {
      result.block(0, i * dim_local_, dim_global_, dim_local_) =
          Da.block(0, i * dim_middle_, dim_global_, dim_middle_) *
          Db.block(0, i * dim_local_, dim_middle_, dim_local_);
    }
    return result;
  }

  Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override {
    Eigen::MatrixXd result;
    if (dim_global_ == dim_middle_ && dim_middle_ == dim_local_) {
      result = a_->JacobianInverseGramian(local);
      auto b = b_->JacobianInverseGramian(local);
      for (int i = 0; i < local.cols(); ++i) {
        result.block(0, dim_middle_ * i, dim_global_, dim_middle_) *=
            b.block(0, dim_local_ * i, dim_middle_, dim_local_);
      }
    } else {
      auto Da = a_->Jacobian(local);
      auto Db = b_->Jacobian(local);
      result.resize(DimGlobal(), local.cols() * dim_local_);

      Eigen::MatrixXd temp(dim_local_, dim_local_);

      for (int i = 0; i < local.cols(); ++i) {
        temp =
            Db.block(0, dim_local_ * i, dim_middle_, dim_local_).transpose() *
            Da.block(0, dim_middle_ * i, dim_global_, dim_middle_).transpose() *
            Da.block(0, dim_middle_ * i, dim_global_, dim_middle_) *
            Db.block(0, dim_local_ * i, dim_middle_, dim_local_);
        result.block(0, dim_local_ * i, dim_global_, dim_local_) =
            temp.colPivHouseholderQr()
                .solve((Da.block(0, dim_middle_ * i, dim_global_, dim_middle_) *
                        Db.block(0, dim_local_ * i, dim_middle_, dim_local_))
                           .transpose())
                .transpose();
      }
    }
    return result;
  }
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override {
    Eigen::VectorXd result;
    if (dim_global_ == dim_middle_ && dim_middle_ == dim_local_) {
      result = a_->IntegrationElement(local).cwiseProduct(
          b_->IntegrationElement(local));
    } else {
      auto Da = a_->Jacobian(local);
      auto Db = b_->Jacobian(local);
      result.resize(local.cols());
      Eigen::MatrixXd temp(dim_local_, dim_local_);
      for (int i = 0; i < local.cols(); ++i) {
        temp =
            Db.block(0, dim_local_ * i, dim_middle_, dim_local_).transpose() *
            Da.block(0, dim_middle_ * i, dim_global_, dim_middle_).transpose() *
            Da.block(0, dim_middle_ * i, dim_global_, dim_middle_) *
            Db.block(0, dim_local_ * i, dim_middle_, dim_local_);
        result(i) = std::sqrt(temp.determinant());
      }
    }
    return result;
  }

  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override {
    return std::make_unique<CompositeGeometry>(a_->SubGeometry(0, 0),
                                               b_->SubGeometry(codim, i));
  }

  CompositeGeometry(std::unique_ptr<Geometry>&& a,
                    std::unique_ptr<Geometry>&& b)
      : a_(std::move(a)),
        b_(std::move(b)),
        dim_local_(b_->DimLocal()),
        dim_global_(a_->DimLocal()),
        dim_middle_(b_->DimGlobal()) {
    LF_ASSERT_MSG(a_->DimLocal() == b_->DimGlobal(),
                  "dimensions are not compatible.");
  }

 private:
  std::unique_ptr<Geometry> a_;
  std::unique_ptr<Geometry> b_;
  dim_t dim_local_, dim_global_, dim_middle_;
};

}  // namespace

std::unique_ptr<Geometry> Compose(std::unique_ptr<Geometry>&& a,
                                  std::unique_ptr<Geometry>&& b) {
  return std::make_unique<CompositeGeometry>(std::move(a), std::move(b));
}
}  // namespace lf::geometry
