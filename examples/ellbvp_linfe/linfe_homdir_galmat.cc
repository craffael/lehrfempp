/**
 * @file linfe_homdir_galmat.cc
 * @brief Write lowest-order Lagrangian FE Galerkin matrix to file in Matrix
 * Market format
 * @author Ralf Hiptmair
 * @date  December 2023
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace LHG {
// abbreviations for types
using size_type = lf::base::size_type;
using glb_idx_t = lf::assemble::glb_idx_t;

// Demo from NumCSE course: Converting a sparse matrix into COO format, see
// https://stackoverflow.com/questions/28685877/convert-an-eigen-matrix-to-triplet-form-c
// http://eigen.tuxfamily.org/dox/group__TutorialSparse.html
template <typename Scalar>
std::vector<Eigen::Triplet<Scalar>> convertToTriplets(
    const Eigen::SparseMatrix<Scalar> &A) {
  // Empty vector of triplets to be grown in the following loop
  std::vector<Eigen::Triplet<Scalar>> triplets{};
  // Loop over row/columns (depending on column/row major format
  for (int k = 0; k < A.outerSize(); ++k) {
    // Loop over inner dimension and obtain triplets corresponding
    // to non-zero entries.
    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A, k); it;
         ++it) {
      // Retrieve triplet data from iterator
      triplets.emplace_back(it.row(), it.col(), it.value());
    }
  }
  return triplets;
}

// Write a sparse matrix in Eigen's internal format to file in Matrix Market
// format https://math.nist.gov/MatrixMarket/formats.html
// clang-format off
  
  /*
%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% This ASCII file represents a sparse MxN matrix with L 
% nonzeros in the following Matrix Market format:
%
% +----------------------------------------------+
% |%%MatrixMarket matrix coordinate real general | <--- header line
% |%                                             | <--+
% |% comments                                    |    |-- 0 or more comment lines
% |%                                             | <--+         
% |    M  N  L                                   | <--- rows, columns, entries
% |    I1  J1  A(I1, J1)                         | <--+
% |    I2  J2  A(I2, J2)                         |    |
% |    I3  J3  A(I3, J3)                         |    |-- L lines
% |        . . .                                 |    |
% |    IL JL  A(IL, JL)                          | <--+
% +----------------------------------------------+   
%
% Indices are 1-based, i.e. A(1,1) is the first element.
%
%=================================================================================
  5  5  8
    1     1   1.000e+00
    2     2   1.050e+01
    3     3   1.500e-02
    1     4   6.000e+00
    4     2   2.505e+02
    4     4  -2.800e+02
    4     5   3.332e+01
    5     5   1.200e+01
   */

// clang-format on
std::ostream &writeMMformat(std::ostream &out,
                            const Eigen::SparseMatrix<double> &A) {
  const unsigned int M = A.rows();
  const unsigned int N = A.cols();
  const std::vector<Eigen::Triplet<double>> A_coo = convertToTriplets(A);
  out << "%%MatrixMarket matrix coordinate real general" << std::endl;
  out << " " << M << " " << N << " " << A_coo.size() << std::endl;
  out << std::setprecision(16);
  for (const Eigen::Triplet<double> trp : A_coo) {
    out << trp.row() << " " << trp.col() << " " << trp.value() << std::endl;
  }
  return out;
}

void writeMMformat(const std::string filename,
                   const Eigen::SparseMatrix<double> &A) {
  std::cout << "Writing " << A.rows() << " x " << A.cols()
            << " sparse matrix to file " << filename << std::endl;
  std::ofstream outf(filename, std::ios::out);
  if (outf.good()) {
    (void)writeMMformat(outf, A);
  } else {
    throw std::runtime_error("Cannot open output file");
  }
}

// Read 2D hybrid mesh from file
// See lecturedemomesh.cc for comments
std::shared_ptr<const lf::mesh::Mesh> readMeshFromFile(std::string filename) {
  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(factory), filename);
  return reader.mesh();
}

// Build Galerkin matrix
Eigen::SparseMatrix<double> assembleGalMat(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Diffusion coefficients: here identity
  auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    return Eigen::Matrix<double, 2, 2>::Identity();
  };
  // Wrap diffusion coefficient into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};

  // Build finite-element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Reference to current mesh
  const lf::mesh::Mesh &mesh{*(fe_space->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space`
  const size_type N_dofs(dofh.NumDofs());
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Helper object for computation of element matrices
  lf::fe::DiffusionElementMatrixProvider<double, decltype(mf_alpha)>
      elmat_builder(fe_space, mf_alpha);
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  // Dummy right-hand side vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Eliminate degrees of freedom on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 2)};
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&bd_flags, &dofh](glb_idx_t gdof_idx) -> std::pair<bool, double> {
        const lf::mesh::Entity &node{dofh.Entity(gdof_idx)};
        return (bd_flags(node) ? std::make_pair(true, 0.0)
                               : std::make_pair(false, 0.0));
      },
      A, phi);
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  return A_crs;
}

}  // namespace LHG

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
int main(int argc, char **argv) {
  std::cout << "LehrFEM++ demo: Reading mesh from file, "
            << "assembling Galerkin matrix, "
            << "writing matrix to file " << std::endl;
  // Name of the .gmsh file containing the mesh
  std::string meshfile;
  // Name of the output file
  std::string outfile;
  // Default choices resides in the same directory as this source file
  const std::filesystem::path here = __FILE__;
  auto meshfilepath = here.parent_path() / "demomesh.msh";

  // Parse command-line arguments using boost,
  // see https://www.boost.org/doc/libs/1_63_0/doc/html/program_options.html
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
    ("help,h", "This help message")
    ("meshfile,f",
     po::value<std::string>(&meshfile)->default_value(meshfilepath.string()),"Input .gmsh file")
    ("outfile,o",
     po::value<std::string>(&outfile)->default_value(std::string("galmat.mm")),"Output .mm file");
  // clang-format on
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);
  if (vm.count("help") > 0) {
    std::cout << "Usage " << argv[0] << " -f <meshfile>.gmsh -o <outfile>.mm"
              << std::endl;
  } else {
    std::cout << "Read mesh from file " << meshfile << std::endl;
    std::cout << "Save Galerkin matrix to file " << outfile << std::endl;

    auto mesh_p = LHG::readMeshFromFile(meshfile);
    const auto A = LHG::assembleGalMat(mesh_p);
    LHG::writeMMformat(outfile, A);
  }

  return 0L;
}
