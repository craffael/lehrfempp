/**
 * @file
 * @brief We build a tensor product mesh on the unit square and use pointwise
 * refinement to refine it.
 * @author Anian Ruoss
 * @date   2018-10-06 16:37:17
 * @copyright MIT License
 */

#include <boost/program_options.hpp>
#include <functional>
#include <iostream>
#include <vector>

#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/refutils.h>
#include "lf/io/io.h"
#include "lf/mesh/utils/utils.h"

using CodimMeshDataSet_t =
    std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<bool>>;

bool PointInTriangle(const Eigen::MatrixXd &tria_coords,
                     const Eigen::Vector2d &point) {
  // calculate barycentric coordinates of point using affine transformation
  // from reference triangle: phi(x_hat) = alpha + beta * x_hat
  Eigen::Vector2d alpha = tria_coords.col(0);
  Eigen::Matrix2d beta;
  beta << tria_coords.col(1) - tria_coords.col(0),
      tria_coords.col(2) - tria_coords.col(0);

  Eigen::Vector2d loc_coords = beta.inverse() * (point - alpha);

  return (0 <= loc_coords(0) && 0 <= loc_coords(1) && loc_coords.sum() <= 1);
}

CodimMeshDataSet_t MarkMesh(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh_ptr,
    const Eigen::MatrixXd &point) {
  // we use CodimMeshDataSet to store whether an edge has to be marked or not
  CodimMeshDataSet_t marked =
      lf::mesh::utils::make_CodimMeshDataSet<bool>(mesh_ptr, 1, false);

  // loop through all cells to check if it contains the point
  for (const lf::mesh::Entity &cell : mesh_ptr->Entities(0)) {
    const lf::geometry::Geometry *geom_ptr = cell.Geometry();
    lf::base::RefEl ref_el = cell.RefEl();
    const Eigen::MatrixXd &ref_el_coords(ref_el.NodeCoords());
    const Eigen::MatrixXd vtx_coords(geom_ptr->Global(ref_el_coords));

    // mark all edges if the point lies inside the cell
    if (ref_el == lf::base::RefEl::kTria()) {
      if (PointInTriangle(vtx_coords, point)) {
        for (const lf::mesh::Entity &edge : cell.SubEntities(1)) {
          marked->operator()(edge) = true;
        }
      }
    } else if (ref_el == lf::base::RefEl::kQuad()) {
      // split quadrilateral into two triangles and check each separately
      Eigen::MatrixXd tria1 = vtx_coords.block(0, 0, 2, 3);
      Eigen::MatrixXd tria2(2, 3);
      tria2 << vtx_coords.col(0), vtx_coords.col(2), vtx_coords.col(3);

      if (PointInTriangle(tria1, point) || PointInTriangle(tria2, point)) {
        for (const lf::mesh::Entity &edge : cell.SubEntities(1)) {
          marked->operator()(edge) = true;
        }
      }
    } else {
      std::cerr << "unknown cell geometry" << std::endl;
    }
  }

  return marked;
}

int main(int argc, char **argv) {
  // define allowed command line arguments:
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  desc.add_options()("help", "Produce this help message")(
      "num_steps", po::value<size_t>()->default_value(5),
      "Number of refinement steps")(
      "pointwise", po::value<bool>()->default_value(true),
      "Whether to use pointwise refinement or not")(
      "point",
      po::value<std::vector<double>>()->multitoken()->default_value(
          std::vector<double>{.5, .5}, ".5, .5"),
      "Point coordinates in unit square (ignored if pointwise=false)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help") != 0u) {
    std::cout << desc << std::endl;
    return 1;
  }

  size_t num_steps = vm["num_steps"].as<size_t>();
  bool pointwise = vm["pointwise"].as<bool>();
  std::vector<double> point_coords = vm["point"].as<std::vector<double>>();
  Eigen::Vector2d point(point_coords.data());

  using size_type = lf::base::size_type;
  using lf::io::TikzOutputCtrl;

  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // build single-cell tensor product mesh on unit square
  lf::mesh::hybrid2d::TPQuadMeshBuilder builder(mesh_factory_ptr);
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0});
  builder.setTopRightCorner(Eigen::Vector2d{1, 1});
  builder.setNoXCells(1);
  builder.setNoYCells(1);
  std::shared_ptr<lf::mesh::Mesh> mesh_ptr = builder.Build();

  // output mesh information
  const lf::mesh::Mesh &mesh = *mesh_ptr;
  lf::mesh::utils::PrintInfo(mesh, std::cout);
  std::cout << std::endl;

  // build mesh hierarchy
  lf::refinement::MeshHierarchy multi_mesh(mesh_ptr, mesh_factory_ptr);

  // mark edges of cells containing point
  auto marker = [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge,
                   CodimMeshDataSet_t mesh_data) -> bool {
    return mesh_data->operator()(edge);
  };

  for (int step = 0; step < num_steps; ++step) {
    // obtain pointer to mesh on finest level
    const size_type n_levels = multi_mesh.NumLevels();
    std::shared_ptr<const lf::mesh::Mesh> mesh_fine =
        multi_mesh.getMesh(n_levels - 1);

    // print number of entities of various co-dimensions
    std::cout << "Mesh on level " << n_levels - 1 << ": "
              << mesh_fine->NumEntities(2) << " nodes, "
              << mesh_fine->NumEntities(1) << " edges, "
              << mesh_fine->NumEntities(0) << " cells," << std::endl;

    lf::io::writeTikZ(
        *mesh_fine,
        std::string("refinement_mesh") + std::to_string(step) + ".txt",
        TikzOutputCtrl::RenderCells | TikzOutputCtrl::CellNumbering |
            TikzOutputCtrl::VerticeNumbering | TikzOutputCtrl::NodeNumbering |
            TikzOutputCtrl::EdgeNumbering);

    lf::io::writeMatplotlib(*mesh_fine, std::string("refinement_mesh") +
                                            std::to_string(step) + ".csv");

    if (pointwise) {
      CodimMeshDataSet_t marked_mesh = MarkMesh(mesh_fine, point);
      multi_mesh.MarkEdges(std::bind(marker, std::placeholders::_1,
                                     std::placeholders::_2, marked_mesh));
      multi_mesh.RefineMarked();
    } else {
      multi_mesh.RefineRegular();
    }
  }

  // generate MATLAB functions describing all levels of mesh hierarchy
  lf::refinement::WriteMatlab(multi_mesh, "pointwise_refinement");

  return 0;
}
