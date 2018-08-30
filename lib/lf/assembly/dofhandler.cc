/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Rudimentary implementation of a general DOF handler interface
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "dofhandler.h"

namespace lf::assemble {
UniformFEDofHandler::UniformFEDofHandler(std::shared_ptr<lf::mesh::Mesh> mesh,
                                         const LocalStaticDOFs &&locdof)
    : mesh_(std::move(mesh)) {
  // TODO: Initialization of index arrays
}

std::vector<gdof_idx_t> UniformFEDofHandler::GetGlobalDofs(
    const lf::mesh::Entity &entity) const {
  return std::vector<gdof_idx_t>{};
}

size_type UniformFEDofHandler::GetNoDofs(const lf::mesh::Entity &entity) const {
  switch (entity.RefEl()) {
    case lf::base::RefEl::kPoint(): {
      return no_dofs_[0];
    }
    case lf::base::RefEl::kSegment(): {
      return no_dofs_[2];
    }
    case lf::base::RefEl::kTria(): {
      return no_dofs_[3];
    }
    case lf::base::RefEl::kQuad(): {
      return no_dofs_[4];
    }
  }
}

}  // namespace lf::assemble
