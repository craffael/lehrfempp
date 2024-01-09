#include "mesh_factory.h"

#include <spdlog/spdlog.h>

#include <iostream>

#include "hybrid2d.h"
#include "lf/base/base.h"

namespace lf::mesh::hybrid2d {

std::shared_ptr<spdlog::logger>& MeshFactory::Logger() {
  static auto logger =
      lf::base::InitLogger("lf::mesh::hybrid2d::MeshFactory::Logger");
  return logger;
}

MeshFactory::size_type MeshFactory::AddPoint(coord_t coord) {
  LF_ASSERT_MSG(coord.rows() == dim_world_,
                "coord has incompatible number of rows.");
  // Create default geometry object for a point from location vector
  hybrid2d::Mesh::GeometryPtr point_geo =
      std::make_unique<geometry::Point>(coord);
  nodes_.emplace_back(std::move(point_geo));
  return nodes_.size() - 1;
}

MeshFactory::size_type MeshFactory::AddPoint(
    std::unique_ptr<geometry::Geometry>&& geometry) {
  // Note: For the sake of purely topological computations meshes without
  // any geometry information may make sense.
  // Moreover, the location of a point can in principle be deduced from
  // geometry information supplied for edges or cells.
  // Hence the next assertion should be removed in the medium run.
  LF_ASSERT_MSG(geometry != nullptr,
                "No creation of a point without a valid geometry");
  LF_ASSERT_MSG(geometry->DimGlobal() == dim_world_,
                "geometry->DimGlobal() != dim_world_");
  LF_ASSERT_MSG(geometry->RefEl() == lf::base::RefEl::kPoint(),
                "Geometry object must belong to a point");
  nodes_.emplace_back(std::move(geometry));
  return nodes_.size() - 1;
}

MeshFactory::size_type MeshFactory::AddEntity(
    base::RefEl ref_el, const nonstd::span<const size_type>& nodes,
    std::unique_ptr<geometry::Geometry>&& geometry) {
  LF_ASSERT_MSG(ref_el.Dimension() > 0,
                "Use AddPoint() to add a node to a mesh.");
  LF_ASSERT_MSG(ref_el.Dimension() <= 2, "ref_el.Dimension > 2");
  if (geometry != nullptr) {
    LF_ASSERT_MSG(geometry->DimGlobal() == dim_world_,
                  "geometry->DimGlobal() != dim_world_");
    LF_ASSERT_MSG(geometry->RefEl() == ref_el, "ref_el != geometry->RefEl()");
  }

  // Add an edge
  if (ref_el == base::RefEl::kSegment()) {
    std::array<size_type, 2> ns{};
    unsigned char count = 0;
    for (const auto& n : nodes) {
      LF_ASSERT_MSG(n < nodes_.size(),
                    "node " << n
                            << " specified in call to AddEntity must be "
                               "inserted with AddPoint() first.");
      LF_ASSERT_MSG(
          count < 2,
          "ref_el = segment, but nodes contains more than 2 node indices");
      ns[count] = n;
      ++count;
    }
    LF_ASSERT_MSG(count == 2,
                  "ref_el = segment but size of nodes was " << count);
    edges_.emplace_back(ns, std::move(geometry));
    return edges_.size() - 1;
  }  // end insertion of an edge

  // otherwise its an element:
  // LF_ASSERT_MSG(geometry, "Geometry is required for elements (codim=0)");
  std::array<size_type, 4> ns{};
  unsigned char count = 0;
  for (const auto& n : nodes) {
    LF_ASSERT_MSG(count < ref_el.NumNodes(),
                  "ref_el = " << ref_el << ", but nodes contains " << count + 1
                              << "node indices");
    LF_ASSERT_MSG(n < nodes_.size(),
                  " Node " << n << " for " << ref_el.ToString()
                           << "  must be inserted with AddNode() first.");
    ns[count] = n;
    ++count;
  }
  LF_ASSERT_MSG(count == ref_el.NumNodes(),
                "ref_el.NumNodes() = " << ref_el.NumNodes()
                                       << ", but argument nodes contained "
                                       << count << " nodes");

  // if we have added a triangle, set the fourth node to -1.
  if (count == 3) {
    ns[3] = -1;
  }
  elements_.emplace_back(ns, std::move(geometry));
  return elements_.size() - 1;
}

std::shared_ptr<mesh::Mesh> MeshFactory::Build() {
  // DIAGNOSTICS
  if (Logger()->should_log(spdlog::level::trace)) {
    std::stringstream ss;
    PrintLists(ss);
    SPDLOG_LOGGER_TRACE(Logger(), ss.str());
  }

  // Obtain points to new mesh object; the actual construction of the
  // mesh is done by the constructor of that object
  mesh::Mesh* mesh_ptr =
      new hybrid2d::Mesh(dim_world_, std::move(nodes_), std::move(edges_),
                         std::move(elements_), check_completeness_);

  // Clear all information supplied to the MeshFactory object
  nodes_ = hybrid2d::Mesh::NodeCoordList{};  // .clear();
  edges_ = hybrid2d::Mesh::EdgeList{};       // .clear();
  elements_ = hybrid2d::Mesh::CellList{};    // clear();

  return std::shared_ptr<mesh::Mesh>(mesh_ptr);
}

// For diagnostic output
void MeshFactory::PrintLists(std::ostream& o) const {
  o << "hybrid2d::MeshFactory: Internal information" << '\n';
  o << nodes_.size() << " nodes:" << '\n';
  for (std::size_t j = 0; j < nodes_.size(); j++) {
    o << "Node " << j << " at ";
    if (nodes_[j] != nullptr) {
      o << (nodes_[j]->Global(Eigen::Matrix<double, 0, 1>())).transpose()
        << '\n';
    } else {
      o << "with unknown location!" << '\n';
    }
  }
  o << edges_.size() << " edges " << '\n';
  for (std::size_t j = 0; j < edges_.size(); j++) {
    o << "Edge " << j << ": " << edges_[j].first[0] << " <-> "
      << edges_[j].first[1];
    if (edges_[j].second) {
      o << " with geometry";
    }
    o << '\n';
  }

  std::cout << elements_.size() << " cells " << '\n';
  for (std::size_t j = 0; j < elements_.size(); j++) {
    o << "Cell " << j << " : ";
    for (int l = 0; l < 4; l++) {
      if (elements_[j].first[l] != static_cast<size_type>(-1)) {
        o << elements_[j].first[l] << " ";
      }
    }
    if (elements_[j].second) {
      o << " with geometry";
    } else {
      o << "[no geometry]";
    }
    o << '\n';
  }
}

}  // namespace lf::mesh::hybrid2d
