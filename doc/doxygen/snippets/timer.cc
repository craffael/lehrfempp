/**
 * @file
 * @brief Show usage of Autotimer
 * @author Raffael Casagrande
 * @date   2020-10-30 05:31:07
 * @copyright MIT License
 */

#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

void AutoTimer() {
  //! [AutoTimer]
  lf::io::GmshReader reader(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2), "mesh.msh");

  // We want to measure how long it takes to assemble a mass matrix on this
  // mesh:
  auto fes =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(reader.mesh());
  lf::uscalfe::ReactionDiffusionElementMatrixProvider emp(
      fes, lf::mesh::utils::MeshFunctionConstant(1.0),
      lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::COOMatrix<double> sm(fes->LocGlobMap().NumDofs(),
                                     fes->LocGlobMap().NumDofs());
  {
    std::cout << "Start assembly" << std::endl;

    // start measuring the runtime from here:
    lf::base::AutoTimer at;
    lf::assemble::AssembleMatrixLocally(0, fes->LocGlobMap(), fes->LocGlobMap(),
                                        emp, sm);
  }  // upon destruction of the AutoTimer `at`, the measured runtime will be
     // reported to std::cout.
  //! [AutoTimer]
}
