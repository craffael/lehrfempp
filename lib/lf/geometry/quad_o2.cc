/**
 * @file
 * @brief Implementation of second-order parametric quadrilaterals
 * @author Anian Ruoss
 * @date   2018-12-27 21:16:17
 * @copyright MIT License
 */

#include "quad_o2.h"
#include "point.h"

namespace lf::geometry {

QuadO2::QuadO2(Eigen::Matrix<double, Eigen::Dynamic, 8> coords)
    : coords_(std::move(coords)),
      alpha_(coords_.rows()),
      beta_(coords_.rows(), 2),
      gamma_(coords_.rows(), 2),
      delta_(coords_.rows()),
      epsilon_(coords_.rows(), 2),
      gamma_x_2_(coords_.rows(), 2),
      epsilon_x_2_(coords_.rows(), 2) {
  /*
   *  3 - 6 - 2                     D - G - C
   *  |       |                     |       |
   *  7       5          ->         H       F
   *  |       |                     |       |
   *  0 - 4 - 1                     A - E - B
   */
  const Eigen::VectorXd& A = coords_.col(0);
  const Eigen::VectorXd& B = coords_.col(1);
  const Eigen::VectorXd& C = coords_.col(2);
  const Eigen::VectorXd& D = coords_.col(3);
  const Eigen::VectorXd& E = coords_.col(4);
  const Eigen::VectorXd& F = coords_.col(5);
  const Eigen::VectorXd& G = coords_.col(6);
  const Eigen::VectorXd& H = coords_.col(7);

  alpha_ << A;
  beta_ << -3. * A - B + 4. * E, -3. * A - D + 4. * H;
  gamma_ << 2. * (A + B) - 4. * E, 2. * (A + D) - 4. * H;
  delta_ << 5. * A - B - 3. * C - D - 4. * (E - F - G + H);
  epsilon_ << -2. * (A + B - C - D) + 4. * (E - G),
      -2. * (A - B - C + D) - 4. * (F - H);

  // coefficients for Jacobian()
  gamma_x_2_ << 2. * gamma_;
  epsilon_x_2_ << 2. * epsilon_;
}

Eigen::MatrixXd QuadO2::Global(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  auto local_squared = local.array().square();
  return ((beta_ * local) + (gamma_ * local_squared.matrix()) +
          (delta_ * local.row(0).cwiseProduct(local.row(1))) +
          (epsilon_ *
           (local_squared * local.colwise().reverse().array()).matrix()))
             .colwise() +
         alpha_;
}

Eigen::MatrixXd QuadO2::Jacobian(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  auto local_squared_reversed = local.array().square().colwise().reverse();
  auto local_colwise_product = local.row(0).cwiseProduct(local.row(1));

  Eigen::MatrixXd tmp_epsilon(gamma_.rows(), 2 * local.cols());
  Eigen::MatrixXd tmp_epsilon_x_2 = epsilon_x_2_.replicate(1, local.cols());
  Eigen::MatrixXd tmp_gamma(gamma_.rows(), 2 * local.cols());

  for (long i = 0; i < local.cols(); ++i) {
    tmp_epsilon.block(0, 2 * i, tmp_epsilon.rows(), 2) =
        epsilon_.rowwise().reverse().array().rowwise() *
        local_squared_reversed.col(i).transpose().array();
    tmp_epsilon_x_2.block(0, 2 * i, tmp_epsilon_x_2.rows(), 2) *=
        local_colwise_product(i);
    tmp_gamma.block(0, 2 * i, tmp_gamma.rows(), 2) =
        gamma_x_2_.array().rowwise() * local.col(i).transpose().array();
  }

  Eigen::MatrixXd local_reversed = local.colwise().reverse();
  local_reversed.resize(1, local.size());

  return beta_.replicate(1, local.cols()) + tmp_gamma + tmp_epsilon +
         tmp_epsilon_x_2 +
         (local_reversed.replicate(delta_.rows(), 1).array().colwise() *
          delta_.array())
             .matrix();
}

Eigen::MatrixXd QuadO2::JacobianInverseGramian(
    const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  Eigen::MatrixXd jac = Jacobian(local);
  Eigen::MatrixXd jacInvGram(jac.rows(), jac.cols());

  if (DimGlobal() == 2) {
    for (long i = 0; i < local.cols(); ++i) {
      auto jacobian = jac.block(0, 2 * i, 2, 2);
      jacInvGram.block(0, 2 * i, 2, 2) << jacobian(1, 1), -jacobian(1, 0),
          -jacobian(0, 1), jacobian(0, 0);
      jacInvGram.block(0, 2 * i, 2, 2) /= jacobian.determinant();
    }
  } else {
    for (long i = 0; i < local.cols(); ++i) {
      auto jacobian = jac.block(0, 2 * i, jac.rows(), 2);

      auto A = jacobian.col(0);
      auto B = jacobian.col(1);
      auto AB = A.dot(B);

      Eigen::MatrixXd tmp(2, 2);
      tmp << B.dot(B), -AB, -AB, A.dot(A);

      jacInvGram.block(0, 2 * i, jac.rows(), 2) =
          jacobian * tmp / tmp.determinant();
    }
  }

  return jacInvGram;
}

Eigen::VectorXd QuadO2::IntegrationElement(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  Eigen::MatrixXd jac = Jacobian(local);
  Eigen::VectorXd intElem(local.cols());

  if (DimGlobal() == 2) {
    for (long i = 0; i < local.cols(); ++i) {
      intElem(i) = std::abs(jac.block(0, 2 * i, 2, 2).determinant());
    }
  } else {
    for (long i = 0; i < local.cols(); ++i) {
      auto jacobian = jac.block(0, 2 * i, jac.rows(), 2);

      auto A = jacobian.col(0);
      auto B = jacobian.col(1);
      auto AB = A.dot(B);

      intElem(i) = std::sqrt(std::abs(A.dot(A) * B.dot(B) - AB * AB));
    }
  }

  return intElem;
}

std::unique_ptr<Geometry> QuadO2::SubGeometry(dim_t codim, dim_t i) const {
  switch (codim) {
    case 0: {
      LF_ASSERT_MSG(i == 0, "i is out of bounds");
      return std::make_unique<QuadO2>(coords_);
    }
    case 1: {
      LF_ASSERT_MSG(0 <= i && i <= 3, "i is out of bounds");
      return std::make_unique<SegmentO2>(
          (Eigen::Matrix<double, Eigen::Dynamic, 3>(DimGlobal(), 3)
               << coords_.col(i),
           coords_.col((i + 1) % 4), coords_.col(i + 4))
              .finished());
    }
    case 2: {
      LF_ASSERT_MSG(0 <= i && i <= 8, "i is out of bounds");
      return std::make_unique<Point>(coords_.col(i));
    }
    default: {
      LF_VERIFY_MSG(false, "codim is out of bounds")
    }
  }
}

std::vector<std::unique_ptr<Geometry>> QuadO2::ChildGeometry(
    const RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kQuad(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  LF_VERIFY_MSG(codim < 3, "Illegal codim " << codim);

  const double hLattice = 1. / static_cast<double>(ref_pat.LatticeConst());
  std::vector<std::unique_ptr<Geometry>> childGeoPtrs = {};

  // get coordinates of childGeometries
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> childPolygons(
      ref_pat.ChildPolygons(codim));

  const int noChildren = childPolygons.size();
  LF_VERIFY_MSG(
      noChildren == ref_pat.NumChildren(codim),
      "NumChildren " << noChildren << " <-> " << ref_pat.NumChildren(codim));

  // create a geometry object for each child
  for (size_t child = 0; child < noChildren; ++child) {
    // codim == 0: a child triangle/quadrilateral is described by a lattice
    // polygon with six/eight vertices
    // codim == 1: a child segment is described by a polygon with three vertices
    // codim == 2: a child point by a single point ("polygon with one corner")
    LF_VERIFY_MSG(
        childPolygons[child].rows() == 2,
        "childPolygons[child].rows() = " << childPolygons[child].rows());
    LF_VERIFY_MSG(
        (childPolygons[child].cols() == (3 - codim)) ||
            (childPolygons[child].cols() == 4),
        "childPolygons[child].cols() = " << childPolygons[child].cols());

    const Eigen::MatrixXd locChildPolygonCoords(
        hLattice * childPolygons[child].cast<double>());

    switch (codim) {
      case 0: {
        if (childPolygons[child].cols() == 3) {
          Eigen::MatrixXd locChildCoords(locChildPolygonCoords.rows(), 6);
          locChildCoords << locChildPolygonCoords,
              (locChildPolygonCoords.col(0) + locChildPolygonCoords.col(1)) /
                  2.,
              (locChildPolygonCoords.col(1) + locChildPolygonCoords.col(2)) /
                  2.,
              (locChildPolygonCoords.col(2) + locChildPolygonCoords.col(0)) /
                  2.;

          childGeoPtrs.push_back(
              std::make_unique<TriaO2>(Global(locChildCoords)));

        } else if (childPolygons[child].cols() == 4) {
          Eigen::MatrixXd locChildCoords(locChildPolygonCoords.rows(), 8);
          locChildCoords << locChildPolygonCoords,
              (locChildPolygonCoords.col(0) + locChildPolygonCoords.col(1)) /
                  2.,
              (locChildPolygonCoords.col(1) + locChildPolygonCoords.col(2)) /
                  2.,
              (locChildPolygonCoords.col(2) + locChildPolygonCoords.col(3)) /
                  2.,
              (locChildPolygonCoords.col(3) + locChildPolygonCoords.col(0)) /
                  2.;

          childGeoPtrs.push_back(
              std::make_unique<QuadO2>(Global(locChildCoords)));

        } else {
          LF_VERIFY_MSG(false, "childPolygons[" << child << "].cols() = "
                                                << childPolygons[child].cols());
        }

        break;
      }
      case 1: {
        Eigen::MatrixXd locChildCoords(locChildPolygonCoords.rows(), 3);
        locChildCoords << locChildPolygonCoords,
            (locChildPolygonCoords.col(0) + locChildPolygonCoords.col(1)) / 2.;

        childGeoPtrs.push_back(
            std::make_unique<SegmentO2>(Global(locChildCoords)));

        break;
      }
      case 2: {
        childGeoPtrs.push_back(
            std::make_unique<Point>(Global(locChildPolygonCoords)));

        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Illegal co-dimension");
      }
    }
  }

  return childGeoPtrs;
}

}  // namespace lf::geometry
