#include <gtest/gtest.h>

#include <mesh_function_interpolation.h>
#include <mesh_function_velocity.h>

#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/uscalfe.h>

#include <vector>

static std::shared_ptr<const lf::refinement::MeshHierarchy> buildMeshes(
    lf::base::size_type levels) {
  // Build a very simple triangular mesh
  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(std::move(factory));
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(1, 1);
  builder.setNumXCells(1);
  builder.setNumYCells(1);
  const auto mesh = builder.Build();
  auto meshes = std::make_shared<lf::refinement::MeshHierarchy>(
      mesh, std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  for (lf::base::size_type lvl = 0; lvl < levels; ++lvl) {
    meshes->RefineRegular();
  }
  return meshes;
}

TEST(projects_ipdg_stokes_post_processing,
     mesh_function_interpolation_concept) {
  using mf_inner =
      projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                   double>;
  using mf = projects::ipdg_stokes::post_processing::MeshFunctionInterpolation<
      mf_inner>;
  ASSERT_TRUE(std::is_copy_constructible_v<mf>);
  ASSERT_TRUE(std::is_move_constructible_v<mf>);
  constexpr bool b = lf::uscalfe::internal::IsMeshFunctionCallable<mf, void>(0);
  ASSERT_TRUE(b);
  ASSERT_TRUE(lf::uscalfe::isMeshFunction<mf>);
}

TEST(projects_ipdg_stokes_post_processing, mesh_function_interpolation) {
  // Initialize the mesh hierarchy and the FESpace on the coarsest mesh
  const auto meshes = buildMeshes(5);
  const auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(
          meshes->getMesh(0));

  // Initialize the basis coefficients for the tests
  const lf::base::size_type base_vertex_count =
      meshes->getMesh(0)->NumEntities(2);
  std::vector<Eigen::VectorXd> tests;
  for (lf::base::size_type vertex = 0; vertex < base_vertex_count; ++vertex) {
    tests.push_back(Eigen::VectorXd::Zero(base_vertex_count));
    tests.back()[vertex] = 1;
  }

  // Run through all coefficient vectors
  for (const Eigen::VectorXd& coeffs : tests) {
    // Initialize a mesh function on the coarsest mesh
    const lf::uscalfe::MeshFunctionFE<double, double> mf_inner(fe_space,
                                                               coeffs);
    // Construct a interpolated mesh function for each refinement hierarchy
    // and compare it with the one on the coarsest mesh
    for (lf::base::size_type lvl = 1; lvl < meshes->NumLevels(); ++lvl) {
      projects::ipdg_stokes::post_processing::MeshFunctionInterpolation
          mf_int(mf_inner, *meshes, 0, lvl);
      const auto mesh = meshes->getMesh(lvl);
      const auto mesh0 = meshes->getMesh(0);
      for (const auto entity : mesh->Entities(0)) {
        const auto geom = entity->Geometry();
        const Eigen::Vector2d center =
            geom->Global(entity->RefEl().NodeCoords()).rowwise().sum() /
            entity->RefEl().NumNodes();
        const double value_fine =
            mf_int(*entity, Eigen::Vector2d::Constant(1. / 3))[0];
        double value_coarse;
        if (1 - center[0] + center[1] <= 1) {
          // The parent is the lower right triangle
          Eigen::Matrix2d pb;
          pb << 1, 0, -1, 1;
          value_coarse =
              mf_inner(*(mesh0->EntityByIndex(0, 1)), pb * center)[0];
        } else {
          // The parent is the upper left triangle
          Eigen::Matrix2d pb;
          pb << 1, -1, 0, 1;
          value_coarse =
              mf_inner(*(mesh0->EntityByIndex(0, 0)), pb * center)[0];
        }
        ASSERT_NEAR(value_fine, value_coarse, 1e-5)
            << "Position = [" << center[0] << ", " << center[1] << "]\n"
            << "Coeffs = [" << coeffs.transpose() << "]\n";
      }
    }
  }
}
