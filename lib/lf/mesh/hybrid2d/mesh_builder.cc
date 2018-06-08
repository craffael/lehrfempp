#include "mesh_builder.h"
#include "entity.h"

namespace lf::mesh::hybrid2d {

MeshBuilder::size_type MeshBuilder::AddPoint(coord_t coord) {
  LF_ASSERT_MSG(!built_, "Build() already called.");
  LF_ASSERT_MSG(coord.rows() == dim_world_,
                "coord has incompatible number of rows.");
  nodes_.emplace_back(std::move(coord));
  return nodes_.size() - 1;
}

MeshBuilder::size_type MeshBuilder::AddElement(
    const base::ForwardRange<const size_type>& nodes,
    std::unique_ptr<geometry::Geometry>&& geometry) {
  LF_ASSERT_MSG(!built_, "Build() already called.");
  LF_ASSERT_MSG(geometry->DimGlobal() == dim_world_,
                "geometry->DimGlobal() != dim_world_");
  LF_ASSERT_MSG(geometry->DimLocal() == 2, "geometry->DimLocal() != 2");

  std::vector<size_type> ns;
  ns.reserve(4);
  for (auto& n : nodes) {
    LF_ASSERT_MSG(n < nodes_.size(),
                  "node " << n
                          << " specified in call to AddElement must be "
                             "inserted with AddNode() first.");
    ns.push_back(n);
  }
  LF_ASSERT_MSG(geometry->RefEl().NumNodes() == ns.size(),
                "mismatch between number of nodes and RefEl of geometry.");

  elements_.emplace_back(std::move(ns), std::move(geometry));
  return elements_.size() - 1;
}

std::unique_ptr<mesh::Mesh> MeshBuilder::Build() {
  built_ = true;
  return std::make_unique<Mesh>(dim_world_, std::move(nodes_),
                                std::move(elements_));
}

}  // namespace lf::mesh::hybrid2d
