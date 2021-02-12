/**
 * @file
 * @brief Implementation of BrepMeshFactoryTransfinite
 * @author Raffael Casagrande
 * @date   2021/02/10
 * @copyright MIT License
 */

#include "brep_mesh_factory_transfinite.h"

namespace lf::brep::geom {

template <int INDEX>
auto BrepMeshFactoryTransfinite::GetAllBrepsThroughPoints(
    nonstd::span<const size_type> node_indices) {
  // check if there is a curve that goes through all nodes:
  using p_t = std::decay_t<decltype(std::get<INDEX>(node_breps_[0])[0])>;

  std::vector<p_t*> curves;
  for (const auto& n : node_indices) {
    if (n >= node_breps_.size()) {
      continue;
    }
    for (auto& c : std::get<INDEX>(node_breps_[n])) {
      curves.push_back(&c);
    }
  }

  std::stable_sort(curves.begin(), curves.end(),
                   [](const auto& a, const auto& b) {
                     return a->first.get() < b->first.get();
                   });

  std::pair<typename p_t::first_type, Eigen::MatrixXd> result;
  result.first = nullptr;

  if (curves.size() < node_indices.size()) {
    // It's impossible that there is a curve that goes through all nodes
    return result;
  }
  auto begin = curves.begin();
  for (auto current = curves.begin(); current != curves.end(); ++current) {
    if (current + 1 == curves.end() ||
        (*(current + 1))->first != (*begin)->first) {
      if (current - begin + 1 >= node_indices.size()) {
        LF_ASSERT_MSG(current - begin + 1 == node_indices.size(),
                      "internal error.");
        // we found a curve that goes through all points
        result.first = (*current)->first;
        result.second.resize(INDEX + 1, current - begin + 1);
        int count = 0;
        for (; begin != current + 1; ++begin) {
          if constexpr (INDEX == 0) {
            result.second(0, count) = (*begin)->second;
          } else {
            result.second.col(count) = (*begin)->second;
          }
          ++count;
        }
        break;
      }
      begin = current + 1;
    }
  }

  return result;
}

BrepMeshFactoryTransfinite::BrepMeshFactoryTransfinite(
    std::unique_ptr<mesh::MeshFactory> mesh_factory,
    std::shared_ptr<interface::BrepModel> brep_model)
    : mesh_factory_(std::move(mesh_factory)), brep_(std::move(brep_model)) {
  LF_ASSERT_MSG(mesh_factory_->DimMesh() == 2,
                "So far only 2D meshes are supported.");
  LF_ASSERT_MSG(
      mesh_factory_->DimWorld() == 2 || mesh_factory_->DimWorld() == 3,
      "DimWorld=2 or 3 expected.");
}

mesh::MeshFactory::size_type BrepMeshFactoryTransfinite::AddEntity(
    base::RefEl ref_el, const nonstd::span<const size_type>& nodes,
    std::unique_ptr<geometry::Geometry>&& geometry) {
  LF_ASSERT_MSG(nodes.size() == ref_el.NumNodes(), "Wrong number of nodes.");
  LF_ASSERT_MSG(ref_el != base::RefEl::kPoint(), "Add points with AddPoint()");

#ifndef NDEBUG
  // Debug Mode: Check that the corner points are correctly mapped:
  auto global = geometry->Global(ref_el.NodeCoords());
  for (int i = 0; i < global.cols(); ++i) {
    LF_ASSERT_MSG(node_coords_.size() > nodes[i],
                  "node index " << i << " is out of bounds.");
    LF_ASSERT_MSG(node_coords_[nodes[i]].rows() > 0,
                  "Error, node #" << i << " was never added.");
    LF_ASSERT_MSG(
        (global.col(i) - node_coords_[nodes[i]]).norm() <
            (1 + node_coords_[nodes[i]].norm()) * 1e-6,
        "Node #"
            << nodes[i] << " has coordinates "
            << node_coords_[nodes[i]].transpose()
            << ", but the geometry passed to AddEntity() maps this node to "
            << global.col(i).transpose());
  }
#endif

  if (ref_el == base::RefEl::kSegment()) {
    // check if there is a curve that goes through both endpoints:
    auto [curve, param] = GetAllBrepsThroughPoints<0>(nodes);
    if (curve != nullptr) {
      Eigen::RowVector2d curve_param(param(0), param(1));
      if (DimWorld() == 2) {
        curve = std::make_shared<BrepCurve2d>(std::move(curve));
      }
      FixCurveParamsIfPeriodic(*curve, curve_param);
      return mesh_factory_->AddEntity(
          base::RefEl::kSegment(), nodes,
          std::make_unique<BrepCurve>(curve, curve_param));
    }
  } else if (ref_el == base::RefEl::kTria()) {
    // check if there is a surface going through all points:
    auto [surface, surf_param] = GetAllBrepsThroughPoints<1>(nodes);
    LF_VERIFY_MSG(surface == nullptr,
                  "The triangle with nodes "
                      << nodes[0] << "," << nodes[1] << "," << nodes[2]
                      << " is part of a surface, but surface triangle is not "
                         "yet supported.");

    // check if there are curves going through any of the two nodes:
    std::array<std::pair<std::shared_ptr<const interface::BrepGeometry>,
                         Eigen::RowVector2d>,
               3>
        curves;
    Eigen::MatrixXd cp0, cp1, cp2;
    std::tie(curves[0].first, cp0) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[0], nodes[1]});

    std::tie(curves[1].first, cp1) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[1], nodes[2]});
    std::tie(curves[2].first, cp2) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[2], nodes[0]});
    if (curves[0].first != nullptr || curves[1].first != nullptr ||
        curves[2].first != nullptr) {
      LF_ASSERT_MSG(
          curves[0].first != curves[1].first ||
              curves[1].first != curves[2].first,
          "All three bounding curves of this triangle are actually the same "
          "curve. This yields a degenerate triangle with volume 0.");

      // at least one of the three segments lies on a curve
      // -> WE use BrepTriaTransfinitePerronnet
      auto node_coords = geometry->Global(ref_el.NodeCoords());

      if (curves[0].first == nullptr) {
        // replace it with a straight segment:
        if (DimWorld() == 2) {
          curves[0].first = std::make_shared<CurveStraightLine<2>>(
              node_coords.col(0), node_coords.col(1));
        } else {
          curves[0].first = std::make_shared<CurveStraightLine<3>>(
              node_coords.col(0), node_coords.col(1));
        }
        curves[0].second = Eigen::RowVector2d(0, 1);
      } else {
        curves[0].second = cp0;
        if (DimWorld() == 2) {
          curves[0].first = std::make_shared<BrepCurve2d>(curves[0].first);
        }
        FixCurveParamsIfPeriodic(*curves[0].first, curves[0].second);
      }
      if (curves[1].first == nullptr) {
        // replace it with a straight segment:
        if (DimWorld() == 2) {
          curves[1].first = std::make_shared<CurveStraightLine<2>>(
              node_coords.col(1), node_coords.col(2));
        } else {
          curves[1].first = std::make_shared<CurveStraightLine<3>>(
              node_coords.col(1), node_coords.col(2));
        }
        curves[1].second = Eigen::RowVector2d(0, 1);
      } else {
        curves[1].second = cp1;
        if (DimWorld() == 2) {
          curves[1].first = std::make_shared<BrepCurve2d>(curves[1].first);
        }
        FixCurveParamsIfPeriodic(*curves[1].first, curves[1].second);
      }
      if (curves[2].first == nullptr) {
        // replace it with a straight segment:
        if (DimWorld() == 2) {
          curves[2].first = std::make_shared<CurveStraightLine<2>>(
              node_coords.col(2), node_coords.col(0));
        } else {
          curves[2].first = std::make_shared<CurveStraightLine<3>>(
              node_coords.col(2), node_coords.col(0));
        }
        curves[2].second = Eigen::RowVector2d(0, 1);
      } else {
        curves[2].second = cp2;
        if (DimWorld() == 2) {
          curves[2].first = std::make_shared<BrepCurve2d>(curves[2].first);
        }
        FixCurveParamsIfPeriodic(*curves[2].first, curves[2].second);
      }

      return mesh_factory_->AddEntity(
          base::RefEl::kTria(), nodes,
          std::make_unique<BrepTriaTransfinitePerronnet>(curves));
    }
  } else if (ref_el == base::RefEl::kQuad()) {
    // check if the quad is part of a surface:
    auto [surface, surface_param] = GetAllBrepsThroughPoints<1>(nodes);

    LF_VERIFY_MSG(
        surface == nullptr,
        "The quadrilateral with nodes "
            << nodes[0] << "," << nodes[1] << "," << nodes[2] << "," << nodes[3]
            << " is part of a surface, but surface quadrilateral is not "
               "yet supported.");

    // check if there is at least one segment on the boundary:
    std::array<std::pair<std::shared_ptr<const interface::BrepGeometry>,
                         Eigen::RowVector2d>,
               4>
        curves;
    std::tie(curves[0].first, curves[0].second) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[0], nodes[1]});
    std::tie(curves[1].first, curves[1].second) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[1], nodes[2]});
    std::tie(curves[2].first, curves[2].second) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[2], nodes[3]});
    std::tie(curves[3].first, curves[3].second) = GetAllBrepsThroughPoints<0>(
        std::array<size_type, 2>{nodes[3], nodes[0]});

    LF_VERIFY_MSG(
        curves[0].first == nullptr && curves[1].first == nullptr &&
            curves[2].first == nullptr && curves[3].first == nullptr,
        "Transfinite interpolation on Quadrilateral is not yet supported.");
  } else {
    LF_VERIFY_MSG(false, "unsupported RefEl");
  }
  return mesh_factory_->AddEntity(ref_el, nodes, std::move(geometry));
}

mesh::MeshFactory::size_type BrepMeshFactoryTransfinite::FindBrepForPoint(
    coord_t coord, size_type index) {
  LF_ASSERT_MSG(coord.rows() == DimWorld(),
                "coord has illegal number of rows.");
  Eigen::Vector3d c;
  if (DimWorld() == 2) {
    c.topRows(2) = coord;
    c.z() = 0.;
  } else {
    c = coord;
  }
  if (index >= node_breps_.size()) {
    node_breps_.resize(index + 1);
  }
  node_breps_[index].first = brep_->FindCurves(c);
  node_breps_[index].second = brep_->FindSurfaces(c);

#ifndef NDEBUG
  // in Debug mode: store node coordinates:
  if (index >= node_coords_.size()) {
    node_coords_.resize(index + 1);
  }
  node_coords_[index] = coord;

  // check that the node coordinates agree:
  for (auto& [curve, param] : node_breps_[index].first) {
    LF_ASSERT_MSG(
        (curve->Global((Eigen::Matrix<double, 1, 1>() << param).finished())
             .topRows(DimWorld()) -
         coord)
                .norm() < (1 + coord.norm()) * 1e-6,
        "FindCurves returned invalid result.");
  }
#endif

  return index;
}

void BrepMeshFactoryTransfinite::FixCurveParamsIfPeriodic(
    const interface::BrepGeometry& curve, Eigen::RowVector2d& param) {
  auto period = curve.Periods()(0);
  if (period == 0) {
    return;
  }
  // calculate base length:
  auto calc = [&](const Eigen::RowVector2d& param) {
    auto jac =
        curve.Jacobian((param(1) - param(0)) * qr_segment_.Points() +
                       param.col(0).replicate(1, qr_segment_.NumPoints()));
    return std::abs(param(1) - param(0)) *
           jac.colwise().norm().dot(qr_segment_.Weights());
  };
  double l0 = calc(param);

  if (param(0) < param(1)) {
    double l1 = calc(param + Eigen::RowVector2d::UnitX() * period);
    if (l1 < l0) {
      param(0) += period;
    }
  } else {
    double l1 = calc(param + Eigen::RowVector2d::UnitY() * period);
    if (l1 < l0) {
      param(1) += period;
    }
  }
}

const quad::QuadRule BrepMeshFactoryTransfinite::qr_segment_ =
    quad::make_QuadRule(base::RefEl::kSegment(), 10);

}  // namespace lf::brep::geom
