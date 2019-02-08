/**
 * @file
 * @brief Implementation of second-order parametric triangles
 * @author Anian Ruoss
 * @date   2018-12-05 22:37:17
 * @copyright MIT License
 */

#include "tria_o2.h"
#include "point.h"

namespace lf::geometry {

TriaO2::TriaO2(Eigen::Matrix<double, Eigen::Dynamic, 6> coords)
    : coords_(std::move(coords)),
      alpha_(coords_.rows()),
      beta_(coords_.rows(), 2),
      gamma_(coords_.rows(), 2),
      delta_(coords_.rows()),
      gamma_x_2_(coords_.rows(), 2) {
  /*
   *  2                                C
   *  | \                              | \
   *  5   4              ->            F   E
   *  |     \                          |     \
   *  0 - 3 - 1                        A - D - B
   */
  const Eigen::VectorXd& A = coords_.col(0);
  const Eigen::VectorXd& B = coords_.col(1);
  const Eigen::VectorXd& C = coords_.col(2);
  const Eigen::VectorXd& D = coords_.col(3);
  const Eigen::VectorXd& E = coords_.col(4);
  const Eigen::VectorXd& F = coords_.col(5);

  alpha_ << A;
  beta_ << 4. * D - 3. * A - B, 4. * F - 3. * A - C;
  gamma_ << 2. * (A + B) - 4. * D, 2. * (A + C) - 4. * F;
  delta_ << 4. * (A + E - D - F);

  // coefficient for Jacobian()
  gamma_x_2_ << 2. * gamma_;
}

Eigen::MatrixXd TriaO2::Global(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  return ((beta_ * local) + (gamma_ * local.array().square().matrix()) +
          (delta_ * local.row(0).cwiseProduct(local.row(1))))
             .colwise() +
         alpha_;
}

Eigen::MatrixXd TriaO2::Jacobian(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  Eigen::MatrixXd tmp(gamma_.rows(), 2 * local.cols());

  for (int i = 0; i < local.cols(); ++i) {
    tmp.block(0, 2 * i, tmp.rows(), 2) =
        gamma_x_2_.array().rowwise() * local.col(i).transpose().array();
  }

  Eigen::MatrixXd local_reversed = local.colwise().reverse();
  local_reversed.resize(1, local.size());

  return beta_.replicate(1, local.cols()) + tmp +
         (local_reversed.replicate(delta_.rows(), 1).array().colwise() *
          delta_.array())
             .matrix();
}

Eigen::MatrixXd TriaO2::JacobianInverseGramian(
    const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  Eigen::MatrixXd jac = Jacobian(local);
  Eigen::MatrixXd jacInvGram(jac.rows(), jac.cols());

  if (DimGlobal() == 2) {
    for (int i = 0; i < local.cols(); ++i) {
      auto jacobian = jac.block(0, 2 * i, 2, 2);
      jacInvGram.block(0, 2 * i, 2, 2) << jacobian(1, 1), -jacobian(1, 0),
          -jacobian(0, 1), jacobian(0, 0);
      jacInvGram.block(0, 2 * i, 2, 2) /= jacobian.determinant();
    }
  } else {
    for (int i = 0; i < local.cols(); ++i) {
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

Eigen::VectorXd TriaO2::IntegrationElement(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  Eigen::MatrixXd jac = Jacobian(local);
  Eigen::VectorXd intElem(local.cols());

  if (DimGlobal() == 2) {
    for (int i = 0; i < local.cols(); ++i) {
      intElem(i) = std::abs(jac.block(0, 2 * i, 2, 2).determinant());
    }
  } else {
    for (int i = 0; i < local.cols(); ++i) {
      auto jacobian = jac.block(0, 2 * i, jac.rows(), 2);

      auto A = jacobian.col(0);
      auto B = jacobian.col(1);
      auto AB = A.dot(B);

      intElem(i) = std::sqrt(std::abs(A.dot(A) * B.dot(B) - AB * AB));
    }
  }

  return intElem;
}

std::unique_ptr<Geometry> TriaO2::SubGeometry(dim_t codim, dim_t i) const {
  switch (codim) {
    case 0: {
      LF_ASSERT_MSG(i == 0, "i is out of bounds");
      return std::make_unique<TriaO2>(coords_);
    }
    case 1: {
      LF_ASSERT_MSG(0 <= i && i <= 2, "i is out of bounds");
      return std::make_unique<SegmentO2>(
          (Eigen::Matrix<double, Eigen::Dynamic, 3>(DimGlobal(), 3)
               << coords_.col(i),
           coords_.col((i + 1) % 3), coords_.col(i + 3))
              .finished());
    }
    case 2: {
      LF_ASSERT_MSG(0 <= i && i <= 5, "i is out of bounds");
      return std::make_unique<Point>(coords_.col(i));
    }
    default: { LF_VERIFY_MSG(false, "codim is out of bounds") }
  }
}

std::vector<std::unique_ptr<Geometry>> TriaO2::ChildGeometry(
    const RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kTria(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  LF_VERIFY_MSG(codim < 3, "Illegal codim " << codim);

  const double hLattice = 1. / static_cast<double>(ref_pat.LatticeConst());
  std::vector<std::unique_ptr<Geometry>> childGeoPtrs = {};

  // get coordinates of childGeometries
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> childPolygons(
      ref_pat.ChildPolygons(codim));

  const int noChildren = childPolygons.size();
  LF_VERIFY_MSG(
      noChildren == ref_pat.noChildren(codim),
      "noChildren " << noChildren << " <-> " << ref_pat.noChildren(codim));

  // create a geometry object for each child
  for (size_t child = 0; child < noChildren; ++child) {
    // codim == 0: a child triangle is described by a lattice polygon with six
    // vertices
    // codim == 1: a child segment is described by a polygon with three vertices
    // codim == 2: a child point by a single point ("polygon with one corner")
    LF_VERIFY_MSG(
        childPolygons[child].rows() == 2,
        "childPolygons[child].rows() = " << childPolygons[child].rows());
    LF_VERIFY_MSG(
        childPolygons[child].cols() == (3 - codim),
        "childPolygons[child].cols() = " << childPolygons[child].cols());

    const Eigen::MatrixXd locChildPolygonCoords(
        hLattice * childPolygons[child].cast<double>());

    switch (codim) {
      case 0: {
        Eigen::MatrixXd locChildCoords(locChildPolygonCoords.rows(), 6);
        locChildCoords << locChildPolygonCoords,
            (locChildPolygonCoords.col(0) + locChildPolygonCoords.col(1)) / 2.,
            (locChildPolygonCoords.col(1) + locChildPolygonCoords.col(2)) / 2.,
            (locChildPolygonCoords.col(2) + locChildPolygonCoords.col(0)) / 2.;

        childGeoPtrs.push_back(
            std::make_unique<TriaO2>(Global(locChildCoords)));

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
      default: { LF_VERIFY_MSG(false, "Illegal co-dimension"); }
    }
  }

  return childGeoPtrs;
}

}  // namespace lf::geometry
