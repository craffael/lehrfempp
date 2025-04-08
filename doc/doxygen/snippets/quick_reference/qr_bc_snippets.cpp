// This file contains Doxygen snippets for the quick reference bc document
// It defines a function to hold all code snippets and includes necessary
// imports.

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

void qr_bc_snippets() {
  {
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
    lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1},
                                            {lf::base::RefEl::kSegment(), 1},
                                            {lf::base::RefEl::kTria(), 0},
                                            {lf::base::RefEl::kQuad(), 1}});
    {
      //! [bc Overview]
      lf::mesh::utils::CodimMeshDataSet<bool> bd_flags =
          lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2);
      //! [bc Overview]
    }

    auto g = [](Eigen::Vector2d x) -> double { return x[0] + x[1]; };
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    auto bd_flags =
        lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1);
    //! [bc One function on the whole boundary]
    auto mf_g = lf::mesh::utils::MeshFunctionGlobal(g);

    std::vector<std::pair<bool, double>> boundary_val =
        lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g);
    //! [bc One function on the whole boundary]

    {
      //! [bc One function on the whole boundary 0]
      auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
        return boundary_val[dof_idx];
      };
      //! [bc One function on the whole boundary 0]
    }

    //! [bc Constant value on the boundary]
    auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
      if (bd_flags(dofh.Entity(dof_idx))) {
        return std::make_pair(true, 1.0);  // value 1.0 on the boundary
      } else {
        return std::make_pair(false, 0.0);  // value irrelevant
      }
    };
    //! [bc Constant value on the boundary]

    //! [bc BC only on part of the boundary]
    lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags_part(mesh_p, false);
    //! [bc BC only on part of the boundary]

    bool condition = true;
    //! [bc BC only on part of the boundary 0]
    for (const auto& edge : fe_space->Mesh()->Entities(1)) {
      if (condition) bd_flags_part(*edge) = true;
    }

    for (const auto& node : fe_space->Mesh()->Entities(2)) {
      //...
    }
    //! [bc BC only on part of the boundary 0]
  }
}
