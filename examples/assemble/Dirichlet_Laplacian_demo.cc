/**
 * @file
 * @brief Implementation of linear Lagrangian finite elements for the Dirichlet
 *        problem for the Laplacian
 * @author Ralf Hiptmair
 * @date   October 2018
 * @copyright MIT License
 */

#include <cmath>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2dp/hybrid2dp.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

/** @brief Assembler model class for linear Lagrangian finite elements
 *
 * The main purpose of this class is to compute the element matrix for
 * the Laplacian on affine triangles or bilinearly mapped quadrilaterals.
 * These element matrices are provided by the `Eval()` method.
 *
 * @note the `Eval()` method will always return a _reference_ to a 4x4 matrix
 * also for triangles. In this case the last row and column must be ignored.
 */
class LinearFELaplaceElementMatrix {
 public:
  using elem_mat_t = Eigen::Matrix<double, 4, 4>;
  using ElemMat = const elem_mat_t &;

  LinearFELaplaceElementMatrix(const lf::mesh::Mesh &mesh) : mesh_(mesh) {}
  bool isActive(const lf::mesh::Entity &) { return true; }
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  elem_mat_t mat_;
  const lf::mesh::Mesh &mesh_;
};

LinearFELaplaceElementMatrix::ElemMat LinearFELaplaceElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  const Eigen::MatrixXd &ref_el_corners(ref_el.NodeCoords());
  const Eigen::MatrixXd vertices{geo_ptr->Global(ref_el_corners)};

  // Computations differ depending on the type of the cell
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      LF_ASSERT_MSG((vertices.cols() == 3) && (vertices.rows() == 2),
                    "Wrong size of vertex matrix!");
      // Obtain area
      Eigen::Matrix<double, 2, 1> refc;
      refc << 1.0 / 3, 1.0 / 3;
      const double area = 0.5 * (geo_ptr->IntegrationElement(refc))[0];

      // Compute the gradients of the barycentric coordinate functions
      // and store them in the columns of a 2x3 matrix grad_bary_coords
      Eigen::Matrix<double, 3, 3> X;  // temporary matrix
      X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
      X.block<3, 2>(0, 1) = vertices.transpose();
      Eigen::Matrix<double, 2, 3> grad_bary_coords{
          X.inverse().block<2, 3>(1, 0)};
      // Since the gradients are constant local integration is easy
      mat_.block<3, 3>(0, 0) =
          area * grad_bary_coords.transpose() * grad_bary_coords;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      LF_ASSERT_MSG((vertices.cols() == 4) && (vertices.rows() == 2),
                    "Wrong size of vertex matrix!");

      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch

  // Return the element matrix
  return mat_;
}

/** @brief Class for computation of local local vector for linear finite
 * elements.
 *
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<double(const Eigen::Vector2d &)>, which supplies
 *         the source function
 *
 * Computation is based on vertex based quadrature
 *
 * @note The element vector returned by the `Eval()` method will always have
 *       length 4 also for triangles
 */

template <typename FUNCTOR>
class LinearFELocalLoadVector {
 public:
  using elem_vec_t = Eigen::Matrix<double, 4, 1>;
  using ElemVec = const elem_vec_t &;

  LinearFELocalLoadVector(const lf::mesh::Mesh &mesh, FUNCTOR f)
      : mesh_(mesh), f_(f) {}
  bool isActive(const lf::mesh::Entity &) { return true; }
  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  elem_vec_t vec_;
  const lf::mesh::Mesh &mesh_;
  FUNCTOR f_;
};

template <typename FUNCTOR>
typename LinearFELocalLoadVector<FUNCTOR>::ElemVec
LinearFELocalLoadVector<FUNCTOR>::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimGlobal() == 2) && (geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  const Eigen::MatrixXd &ref_el_corners(ref_el.NodeCoords());
  const Eigen::MatrixXd vertices{geo_ptr->Global(ref_el_corners)};
  const double area = lf::geometry::Volume(*geo_ptr);

  // get function values in the vertices
  for (int k = 0; k < num_nodes; k++) {
    vec_[k] = area * f_(vertices.col(k)) / num_nodes;
  }
  return vec_;
}

/** @brief Solves Dirichlet problem for the Laplacian
 *
 * @tparam SOLFUNC functor providing boundary values/exact solution
 * @tparam RHSFUNC functor for right hand side source function
 * @param mesh reference to the mesh on which the FE solution is to be computed
 * @param u functor object for solution, also used to sample Dirichlet data
 * @param f object supplying right-hand-side source function
 *
 * @return L2 norm of the error
 */
template <typename SOLFUNC, typename RHSFUNC>
double L2ErrorLinearFEDirichletLaplacian(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, SOLFUNC &&u, RHSFUNC &&f) {
  // Initialize objects for local computations
  LinearFELaplaceElementMatrix loc_mat_laplace(*mesh_p);
  LinearFELocalLoadVector loc_vec_sample(*mesh_p, f);
  // Initialization of index mapping for linear finite elements
  lf::assemble::UniformFEDofHandler loc_glob_map(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});
  // Dimension of finite element space
  const lf::assemble::size_type N_dofs(loc_glob_map.NoDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);
  // Building the Galerkin matrix (trial space = test space)
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, loc_glob_map, loc_mat_laplace);
  // Filling the right-hand-side vector
  auto rhsvec = lf::assemble::AssembleVectorLocally<Eigen::VectorXd>(
      0, loc_glob_map, loc_vec_sample);

  return 0.0;
}

int main(int argc, char **argv) {
  // Pointer to the current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p;

  // Processing command line arguments
  bool verbose = false;
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help,h", "-h -v -f <filename> -s <selection>")
  ("filename,f", po::value<std::string>()->default_value(""),
     "File to load coarse mesh from ")
  ("selector,s", po::value<int>()->default_value(0), "Selection of test mesh")
  ("verbose,v", po::bool_switch(&verbose),"Enable verbose mode");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help") > 0) {
    std::cout << desc << std::endl;
  } else {
    // Retrieve number of degrees of freedom for each entity type from command
    // line arguments
    if (vm.count("filename") > 0) {
      // A filename was specified
      std::string filename{vm["filename"].as<std::string>()};
      if (filename.length() > 0) {
        boost::filesystem::path here = __FILE__;
        auto mesh_file_path = here.parent_path() / filename.c_str();
        auto mesh_factory =
            std::make_unique<lf::mesh::hybrid2dp::MeshFactory>(2);
        lf::io::GmshReader reader(std::move(mesh_factory),
                                  mesh_file_path.string());
        mesh_p = reader.mesh();
      }
    } else {
      if (vm.count("selector") > 0) {
        const int selector = vm["selector"].as<int>();
        if ((selector >= 0) && (selector <= 2)) {
          mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
        }
      }
    }
    if (mesh_p == nullptr) {
      mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
    }
  }
  // At this point a pointer to the mesh is stored in mesh_p

  // Define right hand side
  auto f = [](const Eigen::Vector2d &) { return 0.0; };
  auto u = [](const Eigen::Vector2d &x) {
    return (std::sin(x[0]) * std::sinh(x[1]));
  };
  double L2err = L2ErrorLinearFEDirichletLaplacian(mesh_p,u,f);
}
