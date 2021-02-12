/**
 * @file
 * @brief Demonstrates the use of Transfinite Interpolation elements on a 2d
 * donut and compares it with 2nd order geometry approximations
 * @author Raffael Casagrande
 * @date   2021-02-12 10:32:19
 * @copyright MIT License
 */

#include <lf/brep/geom/geom.h>
#include <lf/brep/occt/occt.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>

#include <boost/program_options.hpp>
#include <filesystem>

using namespace lf;
namespace po = boost::program_options;

void CheckMesh(const mesh::Mesh& mesh, int level) {
  // check the volume:
  double volume = fe::IntegrateMeshFunction(
      mesh, mesh::utils::MeshFunctionConstant(1.0), 20);
  fmt::print("Level {}: Volume: {}, error: {}\n", level, volume,
             (1 - 0.8 * 0.8) * base::kPi - volume);
  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 40);
  //// check if there are elements with very small integration element:
  // for (auto& e : mesh->Entities(0)) {
  //  auto ie = e->Geometry()->IntegrationElement(qr.Points());
  //  if ((ie.array() < 1e-6).any()) {
  //    fmt::print("Level {}: Found an element with determinant {}\n", level,
  //               ie.array().minCoeff());
  //  }
  //}

  // check if there are elements where two edges intersect:
  for (auto& e : mesh.Entities(0)) {
    for (int i = 0; i < 3; ++i) {
      auto edgei = e->SubEntities(1)[i]->Geometry();
      for (int j = i + 1; j < 3; ++j) {
        auto edgej = e->SubEntities(1)[j]->Geometry();
        Eigen::Vector2d x(0.5, 0.5);
        for (int k = 0; k < 100; ++k) {
          Eigen::Matrix2d jac;
          jac.col(0) = edgei->Jacobian(x.row(0));
          jac.col(1) = -edgej->Jacobian(x.row(1));

          Eigen::Vector2d delta = (jac).colPivHouseholderQr().solve(
              edgei->Global(x.row(0)) - edgej->Global(x.row(1)));
          x = x - delta;
          x = x.cwiseMax(0);  // make sure we don't go out of the reference cell
          x = x.cwiseMin(1);
          double dist =
              (edgei->Global(x.row(0)) - edgej->Global(x.row(1))).norm();
          if (delta.norm() < 1e-10 || dist < 1e-10) {
            break;
          }
        }
        if ((std::abs(x(0)) < 1e-6 || std::abs(x(0) - 1) < 1e-6) &&
            (std::abs(x(1)) < 1e-6 || std::abs(x(1) - 1) < 1e-6)) {
          // intersection is at endpoint

        } else {
          double dist =
              (edgei->Global(x.row(0)) - edgej->Global(x.row(1))).norm();
          fmt::print("Min dist: {}\n", dist);
        }
      }
    }
  }
}

int main(int argc, char* argv[]) {
  po::options_description desc;
  desc.add_options()("brep",
                     po::value<std::string>()->default_value("donut2d.brep"),
                     "filename of Brep file")(
      "mesh", po::value<std::string>()->default_value("donut2d.msh"),
      "mesh filename.")("help", "produce this help message");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  // load brep:
  std::string brep_filename = vm["brep"].as<std::string>();
  std::filesystem::path path = __FILE__;
  path = path.parent_path() / brep_filename;
  auto brep_model = std::make_shared<brep::occt::OcctBrepModel>(path.string());

  // Load transfinite mesh:
  std::string msh_filename = vm["mesh"].as<std::string>();
  path = path.parent_path() / msh_filename;
  auto mesh_factory = std::make_unique<brep::geom::BrepMeshFactoryTransfinite>(
      std::make_unique<mesh::hybrid2d::MeshFactory>(2), brep_model);
  auto gmsh_transfinite =
      io::GmshReader(std::move(mesh_factory), path.string());

  // load mesh plainly:
  auto gmsh_plain = io::GmshReader(
      std::make_unique<mesh::hybrid2d::MeshFactory>(2), path.string());

  // refine transfinite mesh:
  refinement::MeshHierarchy mh_transfinite(
      gmsh_transfinite.mesh(),
      std::make_unique<mesh::hybrid2d::MeshFactory>(2));
  int num_level = 2;

  fmt::print("TRANSFINITE INTERPOLATION\n");
  fmt::print("=========================\n");
  for (int level = 0; level < num_level; ++level) {
    if (level > 0) {
      mh_transfinite.RefineRegular();
    }
    CheckMesh(*mh_transfinite.getMesh(level), level);
  }

  // refine plain mesh and check it:
  refinement::MeshHierarchy mh_plain(
      gmsh_plain.mesh(), std::make_unique<mesh::hybrid2d::MeshFactory>(2));
  fmt::print("\n");
  fmt::print("PLAIN MESH               \n");
  fmt::print("=========================\n");
  for (int level = 0; level < num_level; ++level) {
    if (level > 0) {
      mh_plain.RefineRegular();
    }
    CheckMesh(*mh_plain.getMesh(level), level);
  }

  // output mesh to vtk:
  for (auto level = 0; level < mh_transfinite.NumLevels(); ++level) {
    std::string filename = path.stem().string() + ".vtk";

    io::VtkWriter vtk(mh_transfinite.getMesh(level),
                      filename + "_" + std::to_string(level) + ".vtk", 0, 3);
  }
}
