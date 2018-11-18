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
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>
#include "lf/fe/lin_fe.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

static unsigned int dbg_ctrl = 0;
const unsigned int dbg_dofh = 1;
const unsigned int dbg_mesh = 2;
const unsigned int dbg_mat = 4;
const unsigned int dbg_vec = 8;
const unsigned int dbg_bdf = 16;
const unsigned int dbg_elim = 32;
const unsigned int dbg_basic = 64;
const unsigned int dbg_trp = 128;

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
 * No longer in use. Replaced with lf::assemble::fix_flagged_solution_components
 */
void eliminateBoundaryDofs(const std::vector<bool> &tmp_bd_flags,
                           lf::assemble::COOMatrix<double> *A) {
  const lf::base::size_type N = tmp_bd_flags.size();
  LF_ASSERT_MSG((A->cols() == N) && (A->rows() == N),
                "Matrix dimension mismath");

  // Remove rows associated with dofs on the boundary
  auto new_last = std::remove_if(
      A->triplets().begin(), A->triplets().end(),
      [&tmp_bd_flags](lf::assemble::COOMatrix<double>::Triplet &triplet) {
        SWITCHEDSTATEMENT(dbg_ctrl, dbg_elim, if (tmp_bd_flags[triplet.row()]) {
          std::cout << "EBD: removing " << triplet.row() << ',' << triplet.col()
                    << "[" << triplet.value() << "]" << std::endl;
        });
        return tmp_bd_flags[triplet.row()];
      });
  A->triplets().erase(new_last, A->triplets().end());
  // Add unit diagonal entries to rows belonging to dofs on the boundary
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    if (tmp_bd_flags[dofnum]) {
      A->AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

/** @brief Insert sampled values from Dirichlet data into right hand side
 * vector
 *
 */
void insertDirichletDataRHS(const std::vector<bool> &tmp_bd_flags,
                            Eigen::VectorXd *rhs,
                            const Eigen::VectorXd &dirichlet_values) {
  const lf::base::size_type N = tmp_bd_flags.size();
  LF_VERIFY_MSG((rhs->size() == N), "rhs vector size mismatch");
  LF_VERIFY_MSG((dirichlet_values.size() == N), "data vector size mismatch");

  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    if (tmp_bd_flags[dofnum]) {
      // Shape function associated with the boundary
      (*rhs)[dofnum] = dirichlet_values[dofnum];
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
    const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, SOLFUNC &&u,
    RHSFUNC &&f) {
  LF_ASSERT_MSG(mesh_p != nullptr, "Invalid mesh pointer");
  LF_ASSERT_MSG((mesh_p->DimMesh() == 2) && (mesh_p->DimWorld() == 2),
                "For 2D planar meshes only!");
  // Debugging output
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_basic,
      std::cout << "Dirichlet Laplacian: Linear FE L2 error on mesh with "
                << mesh_p->Size(0) << " cells, " << mesh_p->Size(1)
                << " edges, " << mesh_p->Size(2) << " nodes" << std::endl;)
  SWITCHEDSTATEMENT(
      dbg_ctrl, dbg_mesh,
      const int tmp_mesh_ctrl = lf::mesh::hybrid2d::Mesh::output_ctrl_;
      lf::mesh::hybrid2d::Mesh::output_ctrl_ = 100;
      lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
      lf::mesh::hybrid2d::Mesh::output_ctrl_ = tmp_mesh_ctrl);
  // Initialize objects for local computations
  lf::fe::LinearFELaplaceElementMatrix loc_mat_laplace{};
  lf::fe::LinearFELocalLoadVector<double, decltype(f)> loc_vec_sample(f);
  // Initialization of index mapping for linear finite elements
  lf::assemble::UniformFEDofHandler loc_glob_map(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_dofh,
                    std::cout << loc_glob_map << std::endl;);
  // Dimension of finite element space
  const lf::assemble::size_type N_dofs(loc_glob_map.NoDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);
  // Building the Galerkin matrix (trial space = test space)
  // This Galerkin matrix is oblivious of Dirichlet boundary conditions
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, loc_glob_map, loc_mat_laplace);
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_trp, mat.PrintInfo(std::cout));
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_mat,
                    std::cout << "Full " << mat.rows() << 'x' << mat.cols()
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
  std::vector<bool> tmp_bd_flags{flagBoundaryDOFs(loc_glob_map)};
  // >> Old version
  // eliminateBoundaryDofs(tmp_bd_flags, mat);
  // insertDirichletDataRHS(tmp_bd_flags, rhsvec, dirichlet_data);
  // Identify dof indices associated with the boundary
  // >>
  // >> Equivalent new versions
  lf::assemble::fix_flagged_solution_comp_alt<double>(
      [&tmp_bd_flags,
       &dirichlet_data](lf::assemble::gdof_idx_t i) -> std::pair<bool, double> {
        LF_ASSERT_MSG((i < tmp_bd_flags.size()) && (i < dirichlet_data.size()),
                      "Illegal index " << i);
        return std::make_pair(tmp_bd_flags[i], dirichlet_data[i]);
      },
      mat, rhsvec);
  // lf::assemble::fix_flagged_solution_components<double>(
  //   tmp_bd_flags, dirichlet_data, mat, rhsvec);
  // >>
  // Debugging output
  SWITCHEDSTATEMENT(dbg_ctrl, dbg_mat,
                    std::cout << "Reduced " << mat.rows() << 'x' << mat.cols()
                              << " stiffness matrix, " << mat.triplets().size()
                              << " triplets:\n"
                              << mat.makeDense() << std::endl);

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
    const lf::base::size_type num_nodes = cell.RefEl().NumNodes();
    double sum = 0.0;
    for (int k = 0; k < num_nodes; ++k) {
      sum += std::pow(
          sol_vec[cell_dof_idx[k]] - dirichlet_data[cell_dof_idx[k]], 2);
    }
    nodal_err += lf::geometry::Volume(*cell.Geometry()) * (sum / num_nodes);
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
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(std::move(coarse_mesh_p),
                                           mesh_factory_ptr);

  // Perform several steps of regular refinement of the given mesh
  for (int refstep = 0; refstep < reflevels; ++refstep) {
    // Barycentric refinement is the other option
    multi_mesh.RefineRegular(/*lf::refinement::RefPat::rp_barycentric*/);
  }
  // Solve Dirichlet boundary value problem on every level
  lf::assemble::size_type L = multi_mesh.NumLevels();
  std::vector<double> errors(L);
  for (int level = 0; level < L; level++) {
    errors.push_back(
        L2ErrorLinearFEDirichletLaplacian(multi_mesh.getMesh(level), u, f));
  }
  return errors;
}

int main(int argc, const char **argv) {
  // Pointer to the current mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p;

  // Processing command line arguments
  bool verbose = false;
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help,h", "-h -v -f <filename> -s <selection>")
  ("filename,f", "File to load coarse mesh from ")
  ("selector,s", po::value<int>()->default_value(0), "Selection of test mesh")
  ("reflevels,r", po::value<int>()->default_value(2), "Number of refinement levels")
  ("bvpsel,b", po::value<int>()->default_value(0),
   "Selector for Dirichlet data and rhs function")
  ("verbose,v", po::bool_switch(&verbose),"Enable verbose mode");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help") > 0) {
    std::cout << desc << std::endl;
    std::cout << "Internal variables that can be set by name=value args"
              << std::endl;
    lf::base::ListCtrlVars(std::cout);
  } else {
    lf::base::ReadCtrVarsCmdArgs(argc, argv);
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
            std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
        lf::io::GmshReader reader(std::move(mesh_factory),
                                  mesh_file_path.string());
        mesh_p = reader.mesh();
      }
    } else {
      std::cout << "No mesh file supplied, using GenerateHybrid2DTestMesh()"
                << std::endl;
      if (vm.count("selector") > 0) {
        const int selector = vm["selector"].as<int>();
        std::cout << "Using test mesh no " << selector << std::endl;
        if ((selector >= 0) && (selector <= 4)) {
          mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
        }
      }
    }
    // Set number of refinement levels
    unsigned int reflevels = 2;
    if (vm.count("reflevels") > 0) {
      reflevels = vm["reflevels"].as<int>();
    }
    unsigned int bvpsel = 0;
    if (vm.count("bvpsel") > 0) {
      bvpsel = vm["bvpsel"].as<int>();
    }
    if (mesh_p == nullptr) {
      // Default mesh
      std::cout << "Using default mesh; test mesh 0" << std::endl;
      mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
    }
    // At this point a pointer to the mesh is stored in mesh_p
    // Output summary information about the coarsest mesh
    std::cout << "Coarse mesh: " << mesh_p->Size(0) << " cells, "
              << mesh_p->Size(1) << " edges, " << mesh_p->Size(2) << " vertices"
              << std::endl;
    std::cout << reflevels << " refinement levels requested" << std::endl;

    // Problem data provided by function pointers
    std::function<double(const Eigen::Vector2d &)> u, f;

    // Initialize the problem data
    std::cout << "Problem setting " << bvpsel << " selected" << std::endl;
    switch (bvpsel) {
      case 0: {
        // A linear solution, no error, if contained in FE space
        f = [](const Eigen::Vector2d &) { return 0.0; };
        u = [](const Eigen::Vector2d &x) { return (x[0] + 2.0 * x[1]); };
        break;
      }
      case 1: {
        // Quadratic polynomial solution
        f = [](const Eigen::Vector2d &) { return -4.0; };
        u = [](const Eigen::Vector2d &x) {
          return (std::pow(x[0], 2) + std::pow(x[1], 2));
        };
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Illegal problem number");
        break;
      }
    }

    // Set debugging switches
    lf::fe::LinearFELaplaceElementMatrix::dbg_ctrl = 0;
    // LinearFELaplaceElementMatrix::dbg_geo |
    // LinearFELaplaceElementMatrix::dbg_locmat;
    lf::fe::LinearFELocalLoadVector<double, decltype(f)>::dbg_ctrl = 0;
    lf::assemble::DofHandler::output_ctrl_ = 6;
    dbg_ctrl = dbg_basic;  // | dbg_mat | dbg_mesh | dbg_dofh | dbg_trp;
    // lf::assemble::ass_mat_dbg_ctrl = 255;

    // Compute finite element solution and error
    auto L2errs = SolveDirLaplSeqMesh(mesh_p, reflevels, u, f);
    int level = 0;
    for (auto &err : L2errs) {
      std::cout << "L2 error on level " << level << " = " << err << std::endl;
      level++;
    }
  }
  return 0;
}
