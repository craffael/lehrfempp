/**
 * @file
 * @brief Contains the class BrepLagrTria
 * @author Raffael Casagrande
 * @date   2021-01-29 03:45:59
 * @copyright MIT License
 */

#ifndef __33eec0e0a1a440f0a7ad633a16158e14
#define __33eec0e0a1a440f0a7ad633a16158e14

#include <lf/brep/interface/interface.h>
#include <lf/geometry/geometry.h>

#include "brep_lagr_segment.h"
#include "brep_surface_segment.h"

namespace lf::brep::geom {
template <class GEOM_TRIA, class GEOM_SEGMENT>
class BrepLagrTria : public geometry::Geometry {
 public:
  explicit BrepLagrTria(
      GEOM_TRIA tria, const interface::BrepSurface* surface,
      const Eigen::Matrix<double, 2, 3>& surface_coords,
      const std::array<
          std::pair<const interface::BrepCurve*, Eigen::RowVector2d>, 3>&
          curves);

  [[nodiscard]] dim_t DimLocal() const override { return 2; }
  [[nodiscard]] dim_t DimGlobal() const override { return geom_.DimGlobal(); }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override;
  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const geometry::RefinementPattern& ref_pat,
      lf::base::dim_t codim) const override;

 private:
  GEOM_TRIA geom_;
  // Pointer to the surface paramterization, can be null!
  const brep::interface::BrepSurface* surface_;
  // contains the surface coordinates of the tria corners
  Eigen::Matrix<double, 2, 3> surface_coords_;

  // contains the three edge parametrizations (first) and the local coordinates
  // of the endpoints w.r.t. these parametrizations.
  // Note that the vector can be nullptr!
  std::array<std::pair<const interface::BrepCurve*, Eigen::RowVector2d>, 3>
      curves_;
};

template <class GEOM_TRIA, class GEOM_SEGMENT>
BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::BrepLagrTria(
    GEOM_TRIA tria, const interface::BrepSurface* surface,
    const Eigen::Matrix<double, 2, 3>& surface_coords,
    const std::array<std::pair<const interface::BrepCurve*, Eigen::RowVector2d>,
                     3>& curves)
    : geom_(tria),
      surface_(surface),
      surface_coords_(surface_coords),
      curves_(curves) {
  LF_ASSERT_MSG(geom_.RefEl() == base::RefEl::kTria(),
                "BrepLagrTria expects to wrap a tria.");
  // check that the corners of the surface agree with the corners of the
  // geometry
  if (surface != nullptr) {
    LF_ASSERT_MSG(geom_.Global(base::RefEl::kTria().NodeCoords())
                      .isApprox(surface_->GlobalMulti(surface_coords_)),
                  "Corners of Geometry do not agree with surface mapping.");
    LF_ASSERT_MSG(curves_[0].first != nullptr && curves_[1].first != nullptr &&
                      curves_[2].first != nullptr,
                  "If surface geometry is set, also all edge geometries must "
                  "be specified.");
  }
  if (curves_[0].first != nullptr) {
    // check that the corners of the first edge agree with the corners of the
    // geometry
    for (int i = 0; i < 2; ++i) {
      LF_ASSERT_MSG(
          curves_[0]
              .first->GlobalSingle(curves[0].second(i))
              .isApprox(geom_.Global(base::RefEl::kTria().NodeCoords().col(
                  base::RefEl::kTria().SubSubEntity2SubEntity(1, 0, 1, i)))));
    }
  }
  if (curves_[1].first != nullptr) {
    // check that the corners of the first edge agree with the corners of the
    // geometry
    for (int i = 0; i < 2; ++i) {
      LF_ASSERT_MSG(
          curves_[1]
              .first->GlobalSingle(curves[1].second(i))
              .isApprox(geom_.Global(base::RefEl::kTria().NodeCoords().col(
                  base::RefEl::kTria().SubSubEntity2SubEntity(1, 1, 1, i)))));
    }
  }
  if (curves_[2].first != nullptr) {
    // check that the corners of the first edge agree with the corners of the
    // geometry
    for (int i = 0; i < 2; ++i) {
      LF_ASSERT_MSG(
          curves_[2]
              .first->GlobalSingle(curves[2].second(i))
              .isApprox(geom_.Global(base::RefEl::kTria().NodeCoords().col(
                  base::RefEl::kTria().SubSubEntity2SubEntity(1, 2, 1, i)))));
    }
  }
}

template <class GEOM_TRIA, class GEOM_SEGMENT>
Eigen::MatrixXd BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::Global(
    const Eigen::MatrixXd& local) const {
  return geom_.Global(local);
}

template <class GEOM_TRIA, class GEOM_SEGMENT>
Eigen::MatrixXd BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::Jacobian(
    const Eigen::MatrixXd& local) const {
  return geom_.Jacobian(local);
}

template <class GEOM_TRIA, class GEOM_SEGMENT>
Eigen::MatrixXd BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::JacobianInverseGramian(
    const Eigen::MatrixXd& local) const {
  return geom_.JacobianInverseGramian(local);
}

template <class GEOM_TRIA, class GEOM_SEGMENT>
Eigen::VectorXd BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  return geom_.IntegrationElement(local);
}

template <class GEOM_TRIA, class GEOM_SEGMENT>
std::unique_ptr<geometry::Geometry>
BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::SubGeometry(dim_t codim, dim_t i) const {
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "subidx out of range for codim==0.");
    return std::make_unique<BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>>(
        dynamic_cast<const GEOM_TRIA&>(*geom_.SubGeometry(0, 0)),
        surface_coords_, curves_);
  } else if (codim == 1) {
    LF_ASSERT_MSG(i < 3, "subidx out of range for codim==1.");
    return std::make_unique<BrepLagrSegment<GEOM_SEGMENT>>(
        dynamic_cast<const GEOM_SEGMENT&>(*geom_.SubGeometry(1, i)),
        curves_[i].first, curves_[i].second(0), curves_[i].second(1));
  }
  LF_ASSERT_MSG(codim == 2, "codim out of range.");
  LF_ASSERT_MSG(i < 3, "subidx out of range for codim==2");
  return geom_.SubGeometry(2, i);
}

template <class GEOM_TRIA, class GEOM_SEGMENT>
std::vector<std::unique_ptr<geometry::Geometry>>
BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>::ChildGeometry(
    const geometry::RefinementPattern& ref_pat, base::dim_t codim) const {
  auto mapSurfaceLocal = [&](const auto& local) {
    return surface_coords_.col(0) *
               (1 - local.row(0).array() - local.row(1).array()).matrix() +
           surface_coords_.col(1) * local.row(1) +
           surface_coords_.col(2) * local.row(2);
  };
  auto mapSurfaceGlobal = [&](const auto& local) {
    return surface_->GlobalMulti(mapSurfaceLocal(local));
  };
  auto mapEdgeLocal = [&](int index, const auto& local) {
    return curves_[index].second.col(0) * (1 - local.array()).matrix() +
           curves_[index].second.col(1) * local;
  };
  auto mapEdgeGlobal = [&](int index, const auto& local) {
    return curves_[index].first->GlobalMulti(mapEdgeLocal(local));
  };
  auto make11Matrix = [](double value) {
    return (Eigen::Matrix<double, 1, 1>() << value).finished();
  };
  auto lattice_const = ref_pat.LatticeConst();
  double lattice_constd = static_cast<double>(lattice_const);
  LF_ASSERT_MSG(codim >= 0 && codim < 3, "codim out of bounds.");
  std::vector<std::unique_ptr<geometry::Geometry>> result;

  for (auto& e : ref_pat.ChildPolygons(codim)) {
    // map every point:
    Eigen::MatrixXd mappedPoints(DimGlobal(), e.cols());
    for (int i = 0; i < e.cols(); ++i) {
      if (e(0, i) == 0 && e(1, i) == 0) {
        // node 0:
        mappedPoints.col(i) = geom_->Global(Eigen::Vector2d::Zero());
      } else if (e(0, i) == lattice_const && e(1, i) == 0) {
        // node 1:
        mappedPoints.col(i) = geom_->Global(Eigen::Vector2d(1, 0));
      } else if (e(0, i) == 0 && e(1, i) == lattice_const) {
        // node 2:
        mappedPoints.col(i) = geom_->Global(Eigen::Vector2d(0, 1));
      } else if (e(1, i) == 0 && curves_[0].first != nullptr) {
        // edge 0:
        mappedPoints.col(i) =
            mapEdgeGlobal(0, make11Matrix(e(0, i) / lattice_constd));
      } else if (e(0, i) + e(1, i) == lattice_const &&
                 curves_[1].first != nullptr) {
        // edge 1:
        mappedPoints.col(i) = mapEdgeGlobal(1, e(1, i) / lattice_constd);
      } else if (e(0, i) == 0 && curves_[2].first != nullptr) {
        // edge 2:
        mappedPoints.col(i) = mapEdgeGlobal(2, 1. - e(1, i) / lattice_constd);
      } else if (surface_ != nullptr) {
        // in the interior of the triangle:
        mappedPoints.col(i) =
            mapSurfaceGlobal(e.col(i).cast<double>() / lattice_constd);
      } else {
        mappedPoints.col(i) =
            geom_->Global(e.col(i).cast<double>() / lattice_constd);
      }
    }

    if (codim == 2) {
      LF_ASSERT_MSG(e.cols() == 1, "unexpected #cols.");
      result.push_back(std::make_unique<geometry::Point>(mappedPoints));
    } else if (codim == 1) {
      LF_ASSERT_MSG(e.cols() == 2, "unexpected # cols.");
      if (e(1, 0) == 0 && e(1, 1) == 0 && curves_[0].first != nullptr) {
        // edge0:
        result.push_back(std::make_unique<BrepLagrSegment<GEOM_SEGMENT>>(
            GEOM_SEGMENT(std::move(mappedPoints)), curves_[0].first,
            mapEdgeLocal(0, e.row(0).cast<double>() / lattice_constd)));
      } else if (e(0, 0) + e(1, 0) == lattice_const &&
                 e(0, 1) + e(1, 1) == lattice_const &&
                 curves_[1].first != nullptr) {
        // edge 1:
        result.push_back(std::make_unique<BrepLagrSegment<GEOM_SEGMENT>>(
            GEOM_SEGMENT(std::move(mappedPoints))))
      }
    }
  }

  //////////////////////

  if (codim == 2) {
    auto nodes = ref_pat.ChildPolygons(2);
    for (auto& n : nodes) {
      LF_ASSERT_MSG(n.cols() == 1, "unexpected #cols.");
      result.push_back(std::make_unique<geometry::Point>(
          mapSurfaceGlobal(n.cast<double>() / lattice_const)));
    }
  } else if (codim == 1) {
    for (auto& e : ref_pat.ChildPolygons(1)) {
      LF_ASSERT_MSG(e.cols() == 2, "unexpected # cols.");
      // Check if the edge is a sub entity of an existing edge:
      if (e(1, 0) == 0 && e(1, 1) == 0) {
        // sub entity of edge0:
        auto endpoint_local = mapEdgeLocal(0, ) result.push_back(
            std::make_unique<BrepLagrSegment<GEOM_SEGMENT>>(
                GEOM_SEGMENT(mapEdgeGlobal(0, e.row(0)), )));
      }

      auto local =
          ((e.col(0) * (1 - GEOM_SEGMENT::LagrangeNodes().array()).matrix() +
            e.col(1) * GEOM_SEGMENT::LagrangeNodes()) /
           lattice_const)
              .eval();
      result.push_back(
          std::make_unique<BrepLagrSegment<GEOM_SEGMENT>>(GEOM_SEGMENT(
              mapSurfaceGlobal(local),
              new BrepSurfaceSegment(
                  surface_, mapSurfaceLocal(e.cast<double>() / lattice_const)),
              0, 1.)));
    }
  } else if (codim == 0) {
    for (auto& e : ref_pat.ChildPolygons(0)) {
      LF_ASSERT_MSG(e.cols() == 3, "unexpected #cols");
      auto local =
          (((e.col(0) * (1 - GEOM_TRIA::LagrangeNodes().row(0).array() -
                         GEOM_TRIA::LagrangeNode().row(1).array()))
                .matrix() +
            e.col(1) * GEOM_TRIA::LagrangeNodes().row(0) +
            e.col(2) * GEOM_TRIA::LagrangeNodes().row(1)) /
           lattice_const)
              .eval();
      auto local_edge0 =
          ((e.col(0) * (1 - GEOM_SEGMENT::LagrangeNodes().array()).matrix() +
            e.col(1) * GEOM_SEGMENT::LagrangeNodes()) /
           lattice_const)
              .eval();
      auto local_edge1 =
          ((e.col(1) * (1 - GEOM_SEGMENT::LagrangeNodes().array()).matrix() +
            e.col(2) * GEOM_SEGMENT::LagrangeNodes()) /
           lattice_const)
              .eval();
      auto local_edge2 =
          ((e.col(2) * (1 - GEOM_SEGMENT::LagrangeNodes().array()).matrix() +
            e.col(0) * GEOM_SEGMENT::LagrangeNodes()) /
           lattice_const)
              .eval();

      std::array<std::pair<const interface::BrepCurve*, Eigen::RowVector2d>, 3>
          curves;
      curves[0].first =

          result.push_back(std::make_unique <
                           BrepLagrTria<GEOM_TRIA, GEOM_SEGMENT>(
                               GEOM_TRIA(mapSurfaceGlobal(local)), surface_, ))
    }
  }
}
}  // namespace lf::brep::geom

#endif  // __33eec0e0a1a440f0a7ad633a16158e14
