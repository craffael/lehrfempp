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
#include <lf/refinement/refinement.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

/** @brief Local assembler model class for linear Lagrangian finite elements
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
  using ElemMat = const elem_mat_t;

  LinearFELaplaceElementMatrix(const lf::mesh::Mesh &mesh) : mesh_(mesh) {}
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  ElemMat Eval(const lf::mesh::Entity &cell);

 private:
  /// Reference to current mesh
  const lf::mesh::Mesh &mesh_;
  /// quadrature points on reference quadrilateral
  const double kSqrt3 = 1.0 / std::sqrt(3.0);
  const double kZeta_0 = 0.5 - 0.5 * kSqrt3;
  const double kZeta_1 = 0.5 + 0.5 * kSqrt3;
  const std::array<Eigen::Vector2d, 4> kQuadPoints{
      Eigen::Vector2d{kZeta_0, kZeta_0}, Eigen::Vector2d{kZeta_0, kZeta_1},
      Eigen::Vector2d{kZeta_1, kZeta_0}, Eigen::Vector2d{kZeta_1, kZeta_1}};
  // Gradients of refrence shape functions in rows of a matrix
  static Eigen::Matrix<double, 4, 2> DervRefShapFncts(
      const Eigen::Vector2d &xh);
  // Barycenter of triangle
  const Eigen::Matrix<double, 2, 1> kTriabc{
      (Eigen::Matrix<double, 2, 1>() << (1.0 / 3), (1.0 / 3)).finished()};
  // Point in reference square
  const Eigen::Matrix<double, 2, 1> kQuadpt{
      (Eigen::Matrix<double, 2, 1>() << 0.7, 0.3).finished()};

 public:
  static unsigned int dbg_ctrl;
  const unsigned int dbg_det = 1;
  const unsigned int dbg_locmat = 2;
  const unsigned int dbg_J = 4;
  const unsigned int dbg_geo = 8;
};

unsigned int LinearFELaplaceElementMatrix::dbg_ctrl{0};

inline Eigen::Matrix<double, 4, 2>
LinearFELaplaceElementMatrix::DervRefShapFncts(const Eigen::Vector2d &xh) {
  // clang-format off
  return (Eigen::Matrix<double, 4, 2>(4, 2) <<
	  xh[1] - 1, xh[0] - 1,
	  1 - xh[1],  -xh[0]  ,
	  xh[1]    , xh[0]    ,
	  -xh[1]   , 1 - xh[0]
	  ).finished();
  // clang-format on
}

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
  // Matrix storing corner coordinates in its columns
  const Eigen::MatrixXd vertices{geo_ptr->Global(ref_el_corners)};
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_geo,
                    std::cout << ref_el << ", shape = \n"
                              << vertices << std::endl);

  // 4x4 dense matrix for returning result
  elem_mat_t elem_mat = elem_mat_t::Zero(4, 4);

  // Computations differ depending on the type of the cell
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      LF_ASSERT_MSG((vertices.cols() == 3) && (vertices.rows() == 2),
                    "Wrong size of vertex matrix!");
      // Obtain area (not needed, can also be obtained from determinant
      // of auxiliary matrix
      // Eigen::Matrix<double, 2, 1> refc;
      // refc << 1.0 / 3, 1.0 / 3;
      // const double area = 0.5 * (geo_ptr->IntegrationElement(refc))[0];

      // Set up an auxiliary 3x3-matrix with a leading column 1 and
      // the vertex coordinates in its right 3x2 block
      Eigen::Matrix<double, 3, 3> X;  // temporary matrix
      X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
      X.block<3, 2>(0, 1) = vertices.transpose();
      // The determinant of the auxliriary matrix also supplies the determinant
      const double area = 0.5 * X.determinant();
      SWITCHEDSTATEMENT(dbg_ctrl, dbg_det,
                        std::cout
                            << "Area: " << area << " <-> "
                            << 0.5 * (geo_ptr->IntegrationElement(kTriabc))
                            << std::endl);

      // Compute the gradients of the barycentric coordinate functions
      // and store them in the columns of a 2x3 matrix grad_bary_coords
      Eigen::Matrix<double, 2, 3> grad_bary_coords{
          X.inverse().block<2, 3>(1, 0)};
      // Since the gradients are constant local integration is easy
      elem_mat.block<3, 3>(0, 0) =
          area * grad_bary_coords.transpose() * grad_bary_coords;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      LF_ASSERT_MSG((vertices.cols() == 4) && (vertices.rows() == 2),
                    "Wrong size of vertex matrix!");
      // Coefficient matrix for bilinear mapping to actual cell
      // clang-format off
      Eigen::Matrix<double, 2, 4> G(vertices *
                                    (Eigen::Matrix<double, 4, 4>() <<
				     1, -1, -1, 1,
				     0, 1 , 0, -1,
				     0, 0, 0, 1,
				     0, 0, 1, -1)
                                        .finished());
      // clang-format on 
      // Coefficients for the determinant as a linear function in
      // the reference coordinates
      // clang-format off
      Eigen::Vector3d detc{
          (Eigen::Vector3d() <<
	   G(0, 1) * G(1, 2) - G(1, 1) * G(0, 2),
           G(0, 1) * G(1, 3) - G(1, 1) * G(0, 3),
           G(1, 2) * G(0, 3) - G(0, 2) * G(1, 3))
              .finished()};
      // clang-format on
      // Determinant at a point in the reference element
      auto detDPhi = [&detc](const Eigen::Vector2d &xh) -> double {
        return std::abs(detc[0] + detc[1] * xh[0] + detc[2] * xh[1]);
      };
      SWITCHEDSTATEMENT(
          dbg_ctrl, dbg_det,
          std::cout << "Determinant: " << detDPhi(kQuadpt) << " <-> "
                    << (geo_ptr->IntegrationElement(kQuadpt)) << std::endl);

      // Transposed adjunct matrix of Jacobian of transformation
      auto DPhiadj =
          [&G](const Eigen::Vector2d &xh) -> Eigen::Matrix<double, 2, 2> {
        // clang-format off
        return ((Eigen::Matrix<double, 2, 2>() <<
		 G(1, 2) + G(1, 3) * xh[0] , -G(1, 1) - G(1, 3) * xh[1],
		 -G(0, 2) - G(0, 3) * xh[0], G(0, 1) + G(0, 3) * xh[1])
                    .finished());
        // clang-format on
      };
      // Output for debugging
      SWITCHEDSTATEMENT(
          dbg_ctrl, dbg_J,
          std::cout << "InverseTransposedJacobian:\n"
                    << (DPhiadj(kQuadpt) / detDPhi(kQuadpt)) << " <-> "
                    << (geo_ptr->JacobianInverseGramian(kQuadpt)) << std::endl);

      // For a quadrilateral we have to use numerical quadrature based on
      // a 2-2 tensor-product Gauss rule.
      // Sum over quadrature points (weight= 0.5)
      for (int q = 0; q < 4; ++q) {
        // Obtain reference coordinates of current quadrature point
        const Eigen::Vector2d qp(kQuadPoints[q]);
        // Gradients of shape functions in quadrature points stored in
        // the columns of a matrix
        Eigen::Matrix<double, 2, 4> gradmat(DPhiadj(qp) *
                                            DervRefShapFncts(qp).transpose());
        // Update of element matrix by contribution from current quadrature
        // point, whose entries are dot products of gradients
        elem_mat += gradmat.transpose() * gradmat / detDPhi(qp);
      }
      // Scaling with quadrature weights
      elem_mat *= 0.25;
      break;
    }
    default: { LF_ASSERT_MSG(false, "Illegal cell type"); }
  }  // end switch

  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_locmat, lf::base::size_type nv = vertices.cols();
      std::cout << "Element matrix\n"
                << elem_mat.block(0, 0, nv, nv) << std::endl;
      std::cout << "Row sums = " << elem_mat.block(0, 0, nv, nv).colwise().sum()
                << ",\n col sums = "
                << elem_mat.block(0, 0, nv, nv).rowwise().sum().transpose()
                << std::endl);

  // Return the element matrix
  return elem_mat;
}

/** @brief Class for computation of local load vector for linear finite
 * elements.
 *
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<double(const Eigen::Vector2d &)>, which supplies
 *         the source function
 *
 * Computation is based on vertex based quadrature
 *
 * @note The element vector returned by the `Eval()` method will always
 * have length 4 also for triangles
 */
template <typename FUNCTOR>
class LinearFELocalLoadVector {
 public:
  using elem_vec_t = Eigen::Matrix<double, 4, 1>;
  using ElemVec = const elem_vec_t;

  LinearFELocalLoadVector(const lf::mesh::Mesh &mesh, FUNCTOR f)
      : mesh_(mesh), f_(f) {}
  bool isActive(const lf::mesh::Entity &cell) { return true; }
  ElemVec Eval(const lf::mesh::Entity &cell);

 private:
  const lf::mesh::Mesh &mesh_;
  FUNCTOR f_;

 public:
  static unsigned int dbg_ctrl;
  const unsigned int dbg_locvec = 1;
  const unsigned int dbg_geo = 2;
};

template <typename FUNCTOR>
unsigned int LinearFELocalLoadVector<FUNCTOR>::dbg_ctrl = 0;

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
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_geo,
                    std::cout << ref_el << ", shape = \n"
                              << vertices << std::endl);
  const double area = lf::geometry::Volume(*geo_ptr);

  // Vector for returning element vector

  elem_vec_t elem_vec = elem_vec_t::Zero();
  // get function values in the vertices
  for (int k = 0; k < num_nodes; k++) {
    elem_vec[k] = area * f_(vertices.col(k)) / num_nodes;
  }
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_locvec,
      std::cout << "element vector = " << elem_vec.head(num_nodes).transpose()
                << std::endl);
  return elem_vec;
}

static unsigned int dbg_ctrl = 0;
const unsigned int dbg_dofh = 1;
const unsigned int dbg_mesh = 2;
const unsigned int dbg_mat = 4;
const unsigned int dbg_vec = 8;
const unsigned int dbg_bdf = 16;
const unsigned int dbg_elim = 32;
const unsigned int dbg_basic = 64;

/**
 * @brief Build a vector of boundary flags for degrees of freedom
 *
 * @param dofh Local-to-global mapper managing indexing for shape
 * functions
 * @return boolean vector of boundary flags, whose length agrees with the
 *         number of global shape functions managed by `dofh`.
 *
 * Every global shape functions belongs to a unique mesh entity. If that
 * entity is contained in the boundary of the mesh, then the flag array
 * entry with the index of the global shape function is set to `true`,
 * otherwise to `false`.
 */
std::vector<bool> flagBoundaryDOFs(const lf::assemble::DofHandler &dofh) {
  const lf::base::size_type N{dofh.NoDofs()};
  // Flag all entities on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh())};
  // Run through all global shape functions and check whether
  // they are associated with an entity on the boundary. Store
  // this information in a boolean vector
  std::vector<bool> tmp_bd_flags(N, false);
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const lf::mesh::Entity &dof_entity{dofh.Entity(dofnum)};
    tmp_bd_flags[dofnum] = bd_flags(dof_entity);
    SWITCHEDSTATEMENT(
        dbg_ctrl, dbg_bdf,
        std::cout << "FBD: dof " << dofnum << "@ " << dof_entity << " ["
                  << dofh.Mesh()->Index(dof_entity) << "] ";
        if (tmp_bd_flags[dofnum]) { std::cout << "ON BOUNDARY"; } std::cout
        << std::endl);
  }
  return tmp_bd_flags;
}

/** @brief Eliminate degrees of freedom located on the boundary
 *
 */
void eliminateBoundaryDofs(const lf::assemble::DofHandler &dofh,
                           lf::assemble::COOMatrix<double> &A) {
  const lf::base::size_type N{dofh.NoDofs()};
  LF_ASSERT_MSG((A.cols() == N) && (A.rows() == N), "Matrix dimension mismath");
  // Identify dof indices associated with the boundary
  std::vector<bool> tmp_bd_flags{flagBoundaryDOFs(dofh)};
  // Remove rows associated with dofs on the boundary
  lf::assemble::COOMatrix<double>::TripletVec::iterator new_last =
      std::remove_if(
          A.triplets().begin(), A.triplets().end(),
          [&tmp_bd_flags](
              typename lf::assemble::COOMatrix<double>::Triplet &triplet) {
            SWITCHEDSTATEMENT(
                dbg_ctrl, dbg_elim, if (tmp_bd_flags[triplet.row()]) {
                  std::cout << "EBD: removing " << triplet.row() << ','
                            << triplet.col() << "[" << triplet.value() << "]"
                            << std::endl;
                });
            return tmp_bd_flags[triplet.row()];
          });
  A.triplets().erase(new_last, A.triplets().end());
  // Add unit diagonal entries to rows belonging to dofs on the boundary
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    if (tmp_bd_flags[dofnum]) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/** @brief Insert sampled values from Dirichlet data into right hand side
 * vector
 *
 */
void insertDirichletDataRHS(const lf::assemble::DofHandler &dofh,
                            Eigen::VectorXd &rhs,
                            const Eigen::VectorXd &dirichlet_values) {
  const lf::base::size_type N{dofh.NoDofs()};
  LF_VERIFY_MSG((rhs.size() == N), "rhs vector size mismatch");
  LF_VERIFY_MSG((dirichlet_values.size() == N), "data vector size mismatch");
  // Identify dof indices associated with the boundary
  std::vector<bool> tmp_bd_flags{flagBoundaryDOFs(dofh)};
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    if (tmp_bd_flags[dofnum]) {
      // Shape function associated with the boundary
      rhs[dofnum] = dirichlet_values[dofnum];
    }
  }
}

/** @brief Solves Dirichlet problem for the Laplacian
 *
 * @tparam SOLFUNC functor providing boundary values/exact solution
 * @tparam RHSFUNC functor for right hand side source function
 * @param mesh reference to the mesh on which the FE solution is to be
 * computed
 * @param u functor object for solution, also used to sample Dirichlet
 * data
 * @param f object supplying right-hand-side source function
 *
 * @return L2 norm of the nodal error, which is just the L2 norm of the
 *         discretization error approximated by means of the the 2D
 * trapezoidal rule.
 */
template <typename SOLFUNC, typename RHSFUNC>
double L2ErrorLinearFEDirichletLaplacian(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p, SOLFUNC &&u, RHSFUNC &&f) {
  LF_ASSERT_MSG(mesh_p != nullptr, "Invalid mesh pointer");
  LF_ASSERT_MSG((mesh_p->DimMesh() == 2) && (mesh_p->DimWorld() == 2),
                "For 2D planar meshes only!");
  // Debugging output
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_basic,
      std::cout << "Dirichlet Laplacian: Linear FE L2 error on mesh with "
                << mesh_p->Size(0) << " cells, " << mesh_p->Size(1) << " edges, "
                << mesh_p->Size(2) << " nodes" << std::endl;)
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_mesh,
      const int tmp_mesh_ctrl = lf::mesh::hybrid2dp::Mesh::output_ctrl_;
      lf::mesh::hybrid2dp::Mesh::output_ctrl_ = 100;
      lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
      lf::mesh::hybrid2dp::Mesh::output_ctrl_ = tmp_mesh_ctrl);
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
  // This Galerkin matrix is oblivious of Dirichlet boundary conditions
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, loc_glob_map, loc_mat_laplace);
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_mat,
                    std::cout << "Reduced " << mat.rows() << 'x' << mat.cols()
                              << " stiffness matrix, " << mat.triplets().size()
                              << " tripets:\n"
                              << mat.makeDense() << std::endl);

  // Filling the right-hand-side vector
  auto rhsvec = lf::assemble::AssembleVectorLocally<Eigen::VectorXd>(
      0, loc_glob_map, loc_vec_sample);
  // Sample Dirichlet date from the exact solution
  Eigen::VectorXd dirichlet_data(loc_glob_map.NoDofs());
  const Eigen::Matrix<double, 0, 1> ref_coord{
      Eigen::Matrix<double, 0, 1>::Zero(0, 1)};
  for (const lf::mesh::Entity &node : mesh_p->Entities(mesh_p->DimMesh())) {
    LF_ASSERT_MSG(node.RefEl() == lf::base::RefEl::kPoint(),
                  "Wrong topological type for a node");
    const Eigen::Vector2d point = node.Geometry()->Global(ref_coord);
    const lf::assemble::size_type num_int_dof =
        loc_glob_map.NoInteriorDofs(node);
    LF_ASSERT_MSG(num_int_dof == 1, "Node with " << num_int_dof << " dof");
    const lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t> gsf_idx(
        loc_glob_map.InteriorGlobalDofIndices(node));
    const lf::assemble::gdof_idx_t node_dof_idx = gsf_idx[0];
    dirichlet_data[node_dof_idx] = u(point);
  }

  // modify linear system in order to take into account boundary data
  eliminateBoundaryDofs(loc_glob_map, mat);
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_mat,
                    std::cout << "Reduced " << mat.rows() << 'x' << mat.cols()
                              << " stiffness matrix, " << mat.triplets().size()
                              << " triplets:\n"
                              << mat.makeDense() << std::endl);

  insertDirichletDataRHS(loc_glob_map, rhsvec, dirichlet_data);
  // Initialize sparse matrix
  Eigen::SparseMatrix<double> stiffness_matrix(mat.makeSparse());
  // Solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(stiffness_matrix);
  Eigen::VectorXd sol_vec = solver.solve(rhsvec);
  if (solver.info() != Eigen::Success) {
    std::cout << "solver failed!" << std::endl;
  }

  // Compute the norm of nodal error cell by cell
  double nodal_err = 0.0;
  for (const lf::mesh::Entity &cell : mesh_p->Entities(0)) {
    const lf::base::RandomAccessRange<const lf::assemble::gdof_idx_t>
        cell_dof_idx(loc_glob_map.GlobalDofIndices(cell));
    LF_ASSERT_MSG(loc_glob_map.NoLocalDofs(cell) == cell.RefEl().NumNodes(),
                  "Inconsistent node number");
    double sum = 0.0;
    for (int k = 0; k < cell.RefEl().NumNodes(); ++k) {
      sum += std::pow(
          sol_vec[cell_dof_idx[k]] - dirichlet_data[cell_dof_idx[k]], 2);
    }
    nodal_err += lf::geometry::Volume(*cell.Geometry()) * sum;
  }
  return std::sqrt(nodal_err);
}

/** @brief Solves Dirichlet problem for the Laplacian on a sequence of
 *         regularly refined meshes
 * @param coarse_mesh_p pointer to coarsest mesh
 * @param reflevels number of refinements
 * @param u exact solution of the boundary value problem
 * @param f right hand side source function
 * @return L2 norms of discretization errors on each refinement level
 *
 */
template <typename SOLFUNCTOR, typename RHSFUNCTOR>
std::vector<double> SolveDirLaplSeqMesh(
    std::shared_ptr<lf::mesh::Mesh> coarse_mesh_p, unsigned int reflevels,
    SOLFUNCTOR &&u, RHSFUNCTOR &&f) {
  // Prepare for creating a hierarchy of meshes
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(coarse_mesh_p, mesh_factory_ptr);

  // Perform several steps of regular refinement of the given mesh
  for (int refstep = 0; refstep < reflevels; ++refstep) {
    multi_mesh.RefineRegular(/*lf::refinement::RefPat::rp_barycentric*/);
  }
  // Solve Dirichlet boundary value problem on every level
  std::vector<double> errors{};
  for (int level = 0; level < multi_mesh.NumLevels(); level++) {
    errors.push_back(
        L2ErrorLinearFEDirichletLaplacian(multi_mesh.getMesh(level), u, f));
  }
  return errors;
}

int main(int argc, char **argv) {
  // Pointer to the current mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p;

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
  ("reflevels,r", po::value<int>()->default_value(2), "Number of refinement levels")
  ("verbose,v", po::bool_switch(&verbose),"Enable verbose mode");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help") > 0) {
    std::cout << desc << std::endl;
  } else {
    std::cout << "*** Solving Dirichlet problems for the Laplacian ***"
              << std::endl;
    // Retrieve number of degrees of freedom for each entity type from
    // command line arguments
    if (vm.count("filename") > 0) {
      // A filename was specified
      std::string filename{vm["filename"].as<std::string>()};
      if (filename.length() > 0) {
        std::cout << "Reading mesh from file " << filename << std::endl;
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
    // Set number of refinement levels
    unsigned int reflevels = 2;
    if (vm.count("reflevels") > 0) {
      reflevels = vm["reflevels"].as<int>();
    }
    // Default mesh
    if (mesh_p == nullptr) {
      mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
    }
    // At this point a pointer to the mesh is stored in mesh_p
    // Output summary information about the coarsest mesh
    std::cout << "Coarse mesh: " << mesh_p->Size(0) << " cells, "
              << mesh_p->Size(1) << " edges, " << mesh_p->Size(2) << " vertices"
              << std::endl;
    std::cout << reflevels << " refinement levels requested" << std::endl;

    // Define right hand side
    auto f = [](const Eigen::Vector2d &) {
      // return 4.0;
      return 0.0;
    };
    auto u = [](const Eigen::Vector2d &x) {
      // return (std::sin(x[0]) * std::sinh(x[1]));
      // return (std::pow(x[0], 2) + std::pow(x[1], 2));
      return (x[0] + 2.0 * x[1]);
    };

    // Set debugging switches
    LinearFELaplaceElementMatrix::dbg_ctrl = 0;
    LinearFELocalLoadVector<decltype(f)>::dbg_ctrl = 0;
    dbg_ctrl = dbg_basic;

    // Compute finite element solution and error
    auto L2errs = SolveDirLaplSeqMesh(mesh_p,reflevels,u,f);
    int level = 0;
    for (auto & err : L2errs) {
      std::cout << "L2 rrror on level " << level << " = " << err << std::endl;
      level++;
    }
  }
  return 0;
}
