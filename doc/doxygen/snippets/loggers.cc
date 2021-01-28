#include <lf/assemble/assemble.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <iostream>

namespace lf {
void PrintInfoSnippet() {
  //! [PrintInfoMesh]
  std::shared_ptr<mesh::Mesh> mesh;  // initialize it somehow

  // print information about the mesh entities to std::cout :
  mesh::utils::PrintInfo(std::cout, *mesh, 11);
  //! [PrintInfoMesh]
}

void ChangeLogLevel() {
  //! [ChangeLogLevel]
  lf::assemble::AssembleMatrixLogger()->set_level(spdlog::level::debug);
  //! [ChangeLogLevel]

  //! [ChangeLogLevel2]
  lf::mesh::hybrid2d::Mesh::Logger()->set_level(spdlog::level::trace);
  //! [ChangeLogLevel2]

  //! [ChangeLogLevelAll]
  spdlog::apply_all([](const std::shared_ptr<spdlog::logger>& logger) {
    logger->set_level(spdlog::level::debug);
  });
  //! [ChangeLogLevelAll]
}
}  // namespace lf
