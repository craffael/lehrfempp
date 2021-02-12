/**
 * @file
 * @brief A wrapper around another mesh factory, which replaces geometries with
 * their transfinite counterparts.
 * @author Raffael Casagrande
 * @date   2021/02/10
 * @copyright MIT License
 */

#ifndef __1F26F1DC_6338_494B_BFBF_32BF92B08BF8
#define __1F26F1DC_6338_494B_BFBF_32BF92B08BF8

#include <lf/brep/interface/interface.h>
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>

#include "brep_curve.h"
#include "brep_curve_2d.h"
#include "brep_tria_transfinite.h"
#include "curve_straight_segment.h"

namespace lf::brep::geom {

class BrepMeshFactoryTransfinite : public mesh::MeshFactory {
 public:
  BrepMeshFactoryTransfinite(std::unique_ptr<mesh::MeshFactory> mesh_factory,
                             std::shared_ptr<interface::BrepModel> brep_model);

  [[nodiscard]] dim_t DimWorld() const override {
    return mesh_factory_->DimWorld();
  }
  [[nodiscard]] dim_t DimMesh() const override {
    return mesh_factory_->DimMesh();
  }
  size_type AddPoint(coord_t coord) override {
    return FindBrepForPoint(coord, mesh_factory_->AddPoint(coord));
  }
  size_type AddPoint(std::unique_ptr<geometry::Geometry>&& geometry) override {
    Eigen::Matrix<double, 0, 1> zero;
    return FindBrepForPoint(geometry->Global(zero),
                            mesh_factory_->AddPoint(std::move(geometry)));
  }

  size_type AddEntity(base::RefEl ref_el,
                      const nonstd::span<const size_type>& nodes,
                      std::unique_ptr<geometry::Geometry>&& geometry) override;

  [[nodiscard]] std::shared_ptr<mesh::Mesh> Build() override {
    return mesh_factory_->Build();
  }

 private:
  std::unique_ptr<mesh::MeshFactory> mesh_factory_;
  std::shared_ptr<interface::BrepModel> brep_;
  // entry[i].first contains all curves going through point with index i with
  // their curve parameter
  // entry[i].second contains all surfaces going through point with index i with
  // surface parameter
  std::vector<std::pair<
      std::vector<
          std::pair<std::shared_ptr<const interface::BrepGeometry>, double>>,
      std::vector<std::pair<std::shared_ptr<const interface::BrepGeometry>,
                            Eigen::Vector2d>>>>
      node_breps_;

  // used to calculate length of a curve (for fixing orientation in presence of
  // periodic boundary conditions)
  static const quad::QuadRule qr_segment_;

#ifndef NDEBUG
  // In Debug mode: Store additionally the node coordinates to make sure the
  // geometries that are passed via AddEntity have the correct node coordinates.
  std::vector<Eigen::VectorXd> node_coords_;
#endif

  size_type FindBrepForPoint(coord_t coord, size_type index);

  template <int INDEX>
  auto GetAllBrepsThroughPoints(nonstd::span<const size_type> node_indices);

  static void FixCurveParamsIfPeriodic(const interface::BrepGeometry& curve,
                                       Eigen::RowVector2d& param);
};

}  // namespace lf::brep::geom

#endif  // 1F26F1DC_6338_494B_BFBF_32BF92B08BF8
