/**
 * @file
 * @brief Defines the OcctBrepModel class
 * @author Raffael Casagrande
 * @date   2020-11-06 05:18:16
 * @copyright MIT License
 */

#ifndef __c22f41716f574518a166a992467c0172
#define __c22f41716f574518a166a992467c0172

#include <Bnd_OBB.hxx>
#include <Geom_Curve.hxx>
#include <TopoDS_Shape.hxx>

#include "lf/brep/interface/interface.h"

#include <Geom_Surface.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include "lf/mesh/utils/all_codim_mesh_data_set.h"

namespace lf::brep::occt {

class OcctBrepModel : public interface::BrepModel {
 public:
  explicit OcctBrepModel(std::string_view filename);

  [[nodiscard]] std::vector<
      std::pair<std::unique_ptr<interface::BrepGeometry>, Eigen::MatrixXd>>
  FindGeometries(base::dim_t dim, const Eigen::MatrixXd& global) const override;

  [[nodiscard]] base::size_type NumGeometries(base::dim_t dim) const override {
    LF_ASSERT_MSG(dim == 2 || dim == 1,
                  "OcctBRepModel contains only 1d or 2d BrepGeometries");
    if (dim == 1) {
      return edges_.size();
    } else {
      return faces_.size();
    }
  }

 private:
  TopoDS_Shape shape_;
  std::vector<std::pair<Bnd_OBB, TopoDS_Edge>> edges_;
  std::vector<std::pair<Bnd_OBB, TopoDS_Face>> faces_;
};

}  // namespace lf::brep::occt

#endif  // __c22f41716f574518a166a992467c0172
