/**
 * @file
 * @brief Check that the HierarchicScalarFESpace class works as expected
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/quad/gauss_quadrature.h>

#include <cmath>
#include <fstream>
#include <limits>

namespace lf::fe::test {

TEST(fe_space_hierarchic, legendre_orthogonality) {
  // Maximum degrees of Legendre polynomials to test
  const int N = 20;
  const int M = 20;
  // Compute the inner product of the n-th and m-th Legendre polynomial
  // and compare it to the analytic solution
  for (int n = 0; n < N; ++n) {
    for (int m = 0; m < M; ++m) {
      const int degree = n + m;
      const int num_points = degree / 2 + 1;
      const auto [points, weights] = lf::quad::GaussLegendre(num_points);
      // Integrate the product of the n-th and m-th legendre polynomial
      // to check whether the norm is as expected
      double result = 0;
      for (int i = 0; i < num_points; ++i) {
        const double Ln = lf::fe::legendre(n, points[i]);
        const double Lm = lf::fe::legendre(m, points[i]);
        result += weights[i] * Ln * Lm;
      }
      const double expected = n == m ? 1. / (2 * n + 1) : 0.;
      ASSERT_NEAR(result, expected, 1e-8) << "n=" << n << " m=" << m;
    }
  }
}

TEST(fe_space_hierarchic, ilegendre) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned i = 0; i <= N; ++i) {
      const double x = static_cast<double>(i) / N;
      // Compute the legendre polynomial at x
      const double exact = lf::fe::legendre(p, x);
      // Differentiate the integral using central differences
      const double intp = lf::fe::ilegendre(p + 1, x + eps / 2);
      const double intm = lf::fe::ilegendre(p + 1, x - eps / 2);
      const double approx = (intp - intm) / eps;
      // Compare the two values
      if (std::fabs(exact) < 100) {
        ASSERT_NEAR(approx, exact, 1e-5) << "P=" << p << " x=" << x;
      } else {
        ASSERT_NEAR(approx / exact, 1, 1e-5) << "p=" << p << " x=" << x;
      }
    }
  }
}

TEST(fe_space_hierarchic, legendre_dx) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned i = 0; i <= N; ++i) {
      const double x = static_cast<double>(i) / N;
      // Compute the differentiated legendre polynomial at x
      const double exact = lf::fe::legendre_dx(p, x);
      // Differentiate the polynomial using central differences
      const double evalp = lf::fe::legendre(p, x + eps / 2);
      const double evalm = lf::fe::legendre(p, x - eps / 2);
      const double approx = (evalp - evalm) / eps;
      // Compare the two values
      if (std::fabs(exact) < 100) {
        ASSERT_NEAR(approx, exact, 1e-5)
            << "p=" << p << " x=" << x << " exact=" << exact
            << " approx=" << approx;
      } else {
        ASSERT_NEAR(approx / exact, 1, 1e-5)
            << "p=" << p << " x=" << x << " exact=" << exact
            << " approx=" << approx;
      }
    }
  }
}

TEST(fe_space_hierarchic, jacobi_general) {
  const unsigned N = 1000;
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned alpha = 1; alpha < 10; ++alpha) {
      for (unsigned i = 0; i <= N; ++i) {
        const double x = static_cast<double>(i) / N;
        // Compute the jacobi polynomial at x using the accurate scheme
        const double accurate = lf::fe::jacobi(p, alpha, x);
        // Compute the jacobi polynomial at x using the general scheme
        const double general = lf::fe::jacobi(p, alpha, 0, x);
        // Compare the two values
        ASSERT_NEAR(accurate, general, 1e-5)
            << "P=" << p << " alpha=" << alpha << " x=" << x;
      }
    }
  }
}

TEST(fe_space_hierarchic, jacobi_integral) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 20; ++p) {
    for (unsigned alpha = 1; alpha < 10; ++alpha) {
      for (unsigned i = 0; i <= N; ++i) {
        const double x = static_cast<double>(i) / N;
        // Compute the jacobi polynomial at x
        const double exact = lf::fe::jacobi(p, alpha, x);
        // Differentiate the integral using central differences
        const double intp = lf::fe::ijacobi(p + 1, alpha, x + eps / 2);
        const double intm = lf::fe::ijacobi(p + 1, alpha, x - eps / 2);
        const double approx = (intp - intm) / eps;
        // Compare the two values
        if (std::fabs(exact) < 100) {
          ASSERT_NEAR(approx, exact, 1e-5)
              << "P=" << p << " alpha=" << alpha << " x=" << x;
        } else {
          ASSERT_NEAR(approx / exact, 1, 1e-5)
              << "p=" << p << " alpha=" << alpha << " x=" << x;
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, jacobi_derivative) {
  const unsigned N = 1000;
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  for (unsigned p = 0; p < 10; ++p) {
    for (unsigned alpha = 1; alpha < 10; ++alpha) {
      for (unsigned i = 0; i <= N; ++i) {
        const double x = static_cast<double>(i) / N;
        // Compute the differentiated jacobi polynomial at x
        const double exact = lf::fe::jacobi_dx(p, alpha, x);
        // Differentiate the polynomial using central differences
        const double evalp = lf::fe::jacobi(p, alpha, x + eps / 2);
        const double evalm = lf::fe::jacobi(p, alpha, x - eps / 2);
        const double approx = (evalp - evalm) / eps;
        // Compare the two values
        if (std::fabs(exact) < 100) {
          ASSERT_NEAR(approx, exact, 1e-5)
              << "p=" << p << " alpha=" << alpha << " x=" << x
              << " exact=" << exact << " approx=" << approx;
        } else {
          ASSERT_NEAR(approx / exact, 1, 1e-5)
              << "p=" << p << " alpha=" << alpha << " x=" << x
              << " exact=" << exact << " approx=" << approx;
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, continuity) {
  for (int selector = 0; selector <= 8; ++selector) {
    // Get a hybrid test mesh on [0, 3]^2
    const auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    for (unsigned p = 1; p <= 10; ++p) {
      // Generate a HierarchicScalarFESpace with varying polynomial degrees:
      mesh::utils::AllCodimMeshDataSet<unsigned> degrees(mesh);
      int count = 0;
      for (auto& e : mesh->Entities(0)) {
        degrees(*e) = (count++ % 2) + 2;
      }
      for (auto& e : mesh->Entities(1)) {
        degrees(*e) = (count++ % 2) + 2;
      }
      for (auto& e : mesh->Entities(2)) {
        degrees(*e) = 1;
      }
      const auto fe_space =
          std::make_unique<lf::fe::HierarchicScalarFESpace<double>>(mesh,
                                                                    degrees);
      // Test the continuity of each basis function associated with an edge
      for (auto&& cell : mesh->Entities(0)) {
        const auto sfl_cell = fe_space->ShapeFunctionLayout(*cell);
        const auto cell_refel = cell->RefEl();
        const auto cell_num_nodes = cell_refel.NumNodes();
        const auto cell_nodes = cell_refel.NodeCoords();
        const auto edges = cell->SubEntities(1);
        const auto orient = cell->RelativeOrientations();
        int edge_offset = cell_num_nodes;  // where do the next edge dofs start?
        for (int i = 0; i < edges.size(); ++i) {
          Eigen::RowVectorXd edge_eval_coords =
              Eigen::RowVectorXd::LinSpaced(p + 1, 0, 1);

          const Eigen::MatrixXd cell_eval_coords = [&]() {
            Eigen::MatrixXd result(2, p + 1);
            Eigen::Vector2d node0 = cell_nodes.col((i + 0) % cell_num_nodes);
            Eigen::Vector2d node1 = cell_nodes.col((i + 1) % cell_num_nodes);
            if (orient[i] == mesh::Orientation::negative) {
              std::swap(node0, node1);
            }
            result.row(0) = Eigen::RowVectorXd::Constant(p + 1, node0[0]) +
                            (node1[0] - node0[0]) * edge_eval_coords;
            result.row(1) = Eigen::RowVectorXd::Constant(p + 1, node0[1]) +
                            (node1[1] - node0[1]) * edge_eval_coords;
            return result;
          }();
          auto a = cell->Geometry()->Global(cell_eval_coords);
          auto b = edges[i]->Geometry()->Global(edge_eval_coords);

          const auto edge = edges[i];
          const auto sfl_edge = fe_space->ShapeFunctionLayout(*edge);
          ASSERT_TRUE(sfl_edge);
          const auto rsf_edge =
              sfl_edge->EvalReferenceShapeFunctions(edge_eval_coords);
          const auto rsf_cell =
              sfl_cell->EvalReferenceShapeFunctions(cell_eval_coords);
          // Compare each basis function on the edge
          ASSERT_TRUE(sfl_edge->NumRefShapeFunctions(0) ==
                      sfl_cell->NumRefShapeFunctions(1, i))
              << "selector=" << selector << " p=" << p
              << " cell=" << mesh->Index(*cell) << " edge=" << i << std::endl;
          for (int rsf_idx = 0; rsf_idx < sfl_edge->NumRefShapeFunctions(0);
               ++rsf_idx) {
            const Eigen::RowVectorXd rsf_edge_eval = rsf_edge.row(2 + rsf_idx);
            Eigen::RowVectorXd rsf_cell_eval;
            if (orient[i] == mesh::Orientation::positive) {
              rsf_cell_eval = rsf_cell.row(edge_offset + rsf_idx);
            } else {
              rsf_cell_eval =
                  rsf_cell.row(edge_offset + sfl_edge->NumRefShapeFunctions(0) -
                               1 - rsf_idx);
            }

            const double max_diff =
                (rsf_edge_eval - rsf_cell_eval).array().abs().maxCoeff();
            ASSERT_LT(max_diff, 1e-10)
                << "selector=" << selector << " p=" << p
                << " cell=" << mesh->Index(*cell) << " edge=" << i
                << " rsf_idx=" << rsf_idx << "\nrsf_edge_eval=["
                << rsf_edge_eval << "]\nrsf_cell_eval=[" << rsf_cell_eval << "]"
                << std::endl;
          }
          edge_offset += sfl_edge->NumRefShapeFunctions(0);
        }
      }
    }
  }
}

/**
 * @brief Utility function, that generates a 2 element mesh consisting each of a
 * triangle + quadrilateral such that these two elements share all possible
 * edges with all possible orientations. It then asks for a FESpace on such a
 * mesh and makes sure this space is continuous across the shared edge.
 * @tparam SCALAR The scalar type of the ScalarFESpace that is returned by
 * fes_factory.
 * @param fes_factory A functor which is called for every mesh and should return
 * the corresponding ScalarFESpace (as a shared_ptr).
 */
template <class SCALAR>
void CheckContinuity(std::function<std::shared_ptr<fe::ScalarFESpace<SCALAR>>(
                         std::shared_ptr<mesh::Mesh>)>
                         fes_factory) {
  Eigen::Matrix<double, Eigen::Dynamic, 4> quad_nodes(2, 4),
      quad_nodes_shifted(2, 4);
  Eigen::Matrix<double, Eigen::Dynamic, 3> tria_nodes(2, 3),
      tria_nodes_shifted(2, 3);
  std::array<base::size_type, 3> tria_node_indices{{1, 4, 2}};
  // clang-format off
  quad_nodes << 0, 1, 1, 0,
                0, 0, 1, 1;
  tria_nodes << 1,2,1,
                0,0,1;
  // clang-format on
  for (int quad_rot : {0, 1, 2, 3}) {
    for (int tria_rot : {0, 1, 2}) {
      for (bool invert_orient : {false, true}) {
        // construct the mesh:
        lf::mesh::hybrid2d::MeshFactory factory(2);
        factory.AddPoint(Eigen::Vector2d{0., 0.});
        factory.AddPoint(Eigen::Vector2d{1., 0.});
        factory.AddPoint(Eigen::Vector2d{1., 1.});
        factory.AddPoint(Eigen::Vector2d{0., 1.});
        factory.AddPoint(Eigen::Vector2d{2., 0.});

        // quad
        std::array<base::size_type, 4> quad_node_indices_shifted;
        base::size_type quad_index;
        for (int i = 0; i < 4; ++i) {
          quad_nodes_shifted.col(i) = quad_nodes.col((i + quad_rot) % 4);
          quad_node_indices_shifted[i] = (i + quad_rot) % 4;
        }
        auto quad_geom = std::make_unique<geometry::QuadO1>(quad_nodes_shifted);
        quad_index =
            factory.AddEntity(base::RefEl::kQuad(), quad_node_indices_shifted,
                              std::move(quad_geom));

        // tria
        std::array<base::size_type, 3> tria_node_indices_shifted;
        base::size_type tria_index;
        for (int i = 0; i < 3; ++i) {
          tria_nodes_shifted.col(i) = tria_nodes.col((i + tria_rot) % 3);
          tria_node_indices_shifted[i] = tria_node_indices[(i + tria_rot) % 3];
        }
        auto tria_geom = std::make_unique<geometry::TriaO1>(tria_nodes_shifted);
        tria_index =
            factory.AddEntity(base::RefEl::kTria(), tria_node_indices_shifted,
                              std::move(tria_geom));

        // edge
        Eigen::Matrix<double, Eigen::Dynamic, 2> edge_nodes(2, 2);
        base::size_type edge_index;
        if (invert_orient) {
          edge_nodes << 1, 1, 1, 0;
          edge_index = factory.AddEntity(
              base::RefEl::kSegment(), std::array<base::size_type, 2>{{2, 1}},
              std::make_unique<geometry::SegmentO1>(std::move(edge_nodes)));
        } else {
          edge_nodes << 1, 1, 0, 1;
          edge_index = factory.AddEntity(
              base::RefEl::kSegment(), std::array<base::size_type, 2>{{1, 2}},
              std::make_unique<geometry::SegmentO1>(std::move(edge_nodes)));
        }

        auto mesh = factory.Build();
        auto edge = mesh->EntityByIndex(1, edge_index);
        auto tria = mesh->EntityByIndex(0, tria_index);
        auto quad = mesh->EntityByIndex(0, quad_index);
        auto quad_edge_subindex =
            (5 - quad_rot) % 4;  // subindex of the common edge w.r.t. quad
        auto tria_edge_subindex =
            (2 - tria_rot) % 3;  // subindex of the common edge w.r.t. tria
        ASSERT_EQ(tria->SubEntities(1)[tria_edge_subindex], edge);
        ASSERT_EQ(quad->SubEntities(1)[quad_edge_subindex], edge);

        // acquire the ScalarFESpace
        auto fes = fes_factory(mesh);
        auto edge_fe = fes->ShapeFunctionLayout(*edge);
        auto tria_fe = fes->ShapeFunctionLayout(*tria);
        auto quad_fe = fes->ShapeFunctionLayout(*quad);

        // transform edge local coordinates to quad/tria local coordinates:
        Eigen::RowVectorXd edge_local =
            Eigen::RowVectorXd::LinSpaced(edge_fe->Degree() + 1, 0, 1);
        Eigen::Matrix<double, 2, Eigen::Dynamic> tria_local(2,
                                                            edge_local.cols()),
            quad_local(2, edge_local.cols());
        const auto& tria_local_node_coords = base::RefEl::kTria().NodeCoords();
        const auto& quad_local_node_coords = base::RefEl::kQuad().NodeCoords();
        base::sub_idx_t tria_subsubindex0 =
            base::RefEl::kTria().SubSubEntity2SubEntity(1, tria_edge_subindex,
                                                        1, 0);
        base::sub_idx_t tria_subsubindex1 =
            base::RefEl::kTria().SubSubEntity2SubEntity(1, tria_edge_subindex,
                                                        1, 1);
        base::sub_idx_t quad_subsubindex0 =
            base::RefEl::kQuad().SubSubEntity2SubEntity(1, quad_edge_subindex,
                                                        1, 0);
        base::sub_idx_t quad_subsubindex1 =
            base::RefEl::kQuad().SubSubEntity2SubEntity(1, quad_edge_subindex,
                                                        1, 1);
        if (tria->RelativeOrientations()[tria_edge_subindex] ==
            mesh::Orientation::positive) {
          tria_local =
              tria_local_node_coords.col(tria_subsubindex0) *
                  (1 - edge_local.array()).matrix() +
              tria_local_node_coords.col(tria_subsubindex1) * edge_local;
        } else {
          tria_local =
              tria_local_node_coords.col(tria_subsubindex1) *
                  (1 - edge_local.array()).matrix() +
              tria_local_node_coords.col(tria_subsubindex0) * edge_local;
        }
        if (quad->RelativeOrientations()[quad_edge_subindex] ==
            mesh::Orientation::positive) {
          quad_local =
              quad_local_node_coords.col(quad_subsubindex0) *
                  (1 - edge_local.array()).matrix() +
              quad_local_node_coords.col(quad_subsubindex1) * edge_local;
        } else {
          quad_local =
              quad_local_node_coords.col(quad_subsubindex1) *
                  (1 - edge_local.array()).matrix() +
              quad_local_node_coords.col(quad_subsubindex0) * edge_local;
        }
        ASSERT_TRUE(quad->Geometry()
                        ->Global(quad_local)
                        .isApprox(edge->Geometry()->Global(edge_local)));
        ASSERT_TRUE(tria->Geometry()
                        ->Global(tria_local)
                        .isApprox(edge->Geometry()->Global(edge_local)));

        // go through all shape functions and make sure they are continuous
        // across the edge:
        Eigen::VectorXd coeff_vector(fes->LocGlobMap().NumDofs());
        for (base::size_type dof = 0; dof < fes->LocGlobMap().NumDofs();
             ++dof) {
          coeff_vector.setZero();
          coeff_vector(dof) = 1;

          auto mf_fe = fe::MeshFunctionFE(fes, coeff_vector);
          auto edge_values = mf_fe(*edge, edge_local);
          auto tria_values = mf_fe(*tria, tria_local);
          auto quad_values = mf_fe(*quad, quad_local);

          ASSERT_EQ(edge_values.size(), edge_local.cols());
          ASSERT_EQ(quad_values.size(), edge_local.cols());
          ASSERT_EQ(tria_values.size(), edge_local.cols());

          for (int i = 0; i < edge_values.size(); ++i) {
            ASSERT_NEAR(edge_values[i], quad_values[i], 1e-7);
            ASSERT_NEAR(edge_values[i], tria_values[i], 1e-7);
          }
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, continuity2) {
  // check continuity with uniformly distributed degrees:
  CheckContinuity<double>([](std::shared_ptr<mesh::Mesh> mesh) {
    return std::make_shared<HierarchicScalarFESpace<double>>(mesh, 10);
  });

  // when degrees vary from entity to entity:
  CheckContinuity<double>([](std::shared_ptr<mesh::Mesh> mesh) {
    auto degrees =
        std::make_unique<mesh::utils::AllCodimMeshDataSet<unsigned>>(mesh);
    int count = 0;
    for (auto& e : mesh->Entities(0)) {
      (*degrees)(*e) = (count++ % 2) + 2;
    }
    for (auto& e : mesh->Entities(1)) {
      (*degrees)(*e) = (count++ % 2) + 2;
    }
    for (auto& e : mesh->Entities(2)) {
      (*degrees)(*e) = 1;
    }
    return std::make_shared<lf::fe::HierarchicScalarFESpace<double>>(
        mesh, [degrees{std::move(degrees)}](const mesh::Entity& e) {
          return (*degrees)(e);
        });
  });
}

TEST(fe_space_hierarchic, segment_dual) {
  // Test the nodal values to dofs function of all segments
  // with degree up to 20
  quad::QuadRuleCache qr_cache;
  for (unsigned p = 1; p <= 20; ++p) {
    const lf::fe::FeHierarchicSegment<double> sfl(p, qr_cache);
    // Get the evaluation nodes for the segment
    const auto eval_nodes = sfl.EvaluationNodes();
    // Evaluate the basis functions at the evaluation nodes
    const Eigen::MatrixXd basis = sfl.EvalReferenceShapeFunctions(eval_nodes);
    // Interpolate all basis functions
    for (unsigned i = 0; i < basis.rows(); ++i) {
      // Compute the basis function coefficients
      const auto dofs = sfl.NodalValuesToDofs(basis.row(i));
      // Only the i-th entry should be 1, all others 0
      for (long j = 0; j < dofs.size(); ++j) {
        ASSERT_NEAR(dofs[j], j == i ? 1 : 0, 1e-8)
            << "p=" << p << " i=" << i << " j=" << j << "\ndofs=[" << dofs
            << "]";
      }
    }
  }
}

TEST(fe_space_hierarchic, tria_dual) {
  // Test the nodal values to dofs function of all segments
  // with degree up to 6
  unsigned max_p = 4;
#ifdef NDEBUG
  max_p = 6;
#endif
  quad::QuadRuleCache qr_cache;
  for (unsigned p_interior = 1; p_interior <= max_p; ++p_interior) {
    for (unsigned p_edge0 = 1; p_edge0 <= max_p; ++p_edge0) {
      for (unsigned p_edge1 = 1; p_edge1 <= max_p; ++p_edge1) {
        for (unsigned p_edge2 = 1; p_edge2 <= max_p; ++p_edge2) {
          // Test for all possible combinations of orientations
          for (const auto o0 : {lf::mesh::Orientation::positive,
                                lf::mesh::Orientation::negative}) {
            for (const auto o1 : {lf::mesh::Orientation::positive,
                                  lf::mesh::Orientation::negative}) {
              for (const auto o2 : {lf::mesh::Orientation::positive,
                                    lf::mesh::Orientation::negative}) {
                const lf::mesh::Orientation orientations[] = {o0, o1, o2};
                std::array<unsigned, 3> p_edges = {{p_edge0, p_edge1, p_edge2}};
                const lf::fe::FeHierarchicTria<double> sfl(
                    p_interior, p_edges, qr_cache, orientations);
                // Get the evaluation nodes for the quad
                const auto eval_nodes = sfl.EvaluationNodes();
                // Evaluate the basis functions at the evaluation nodes
                const Eigen::MatrixXd basis =
                    sfl.EvalReferenceShapeFunctions(eval_nodes);
                // Interpolate all basis functions
                for (unsigned i = 0; i < basis.rows(); ++i) {
                  // Compute the basis function coefficients
                  const auto dofs = sfl.NodalValuesToDofs(basis.row(i));
                  // Only the i-th entry should be 1, all others 0
                  for (long j = 0; j < dofs.size(); ++j) {
                    ASSERT_NEAR(dofs[j], j == i ? 1 : 0, 1e-8)
                        << "o=" << to_char(o0) << to_char(o1) << to_char(o2)
                        << " p=" << p_interior << "," << p_edge0 << ","
                        << p_edge1 << "," << p_edge2 << " i=" << i << " j=" << j
                        << "\ndofs=[" << dofs << "]";
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, quad_dual) {
  // Test the nodal values to dofs function of all segments
  // with degree up to 6
  unsigned max_p = 4;
#ifdef NDEBUG
  max_p = 6;
#endif
  quad::QuadRuleCache qr_cache;
  for (unsigned p_interior = 1; p_interior <= max_p; ++p_interior) {
    for (unsigned p_edge0 = 1; p_edge0 <= max_p; ++p_edge0) {
      for (unsigned p_edge1 = 1; p_edge1 <= max_p; ++p_edge1) {
        for (unsigned p_edge2 = 1; p_edge2 <= max_p; ++p_edge2) {
          for (unsigned p_edge3 = 1; p_edge3 <= max_p; ++p_edge3) {
            // Test for all possible combinations of orientations
            for (const auto o0 : {lf::mesh::Orientation::positive,
                                  lf::mesh::Orientation::negative}) {
              for (const auto o1 : {lf::mesh::Orientation::positive,
                                    lf::mesh::Orientation::negative}) {
                for (const auto o2 : {lf::mesh::Orientation::positive,
                                      lf::mesh::Orientation::negative}) {
                  for (const auto o3 : {lf::mesh::Orientation::positive,
                                        lf::mesh::Orientation::negative}) {
                    const lf::mesh::Orientation orientations[] = {o0, o1, o2,
                                                                  o3};
                    std::array<unsigned, 4> edge_degrees{
                        {p_edge0, p_edge1, p_edge2, p_edge3}};
                    const lf::fe::FeHierarchicQuad<double> sfl(
                        p_interior, edge_degrees, qr_cache, orientations);
                    // Get the evaluation nodes for the quad
                    const auto eval_nodes = sfl.EvaluationNodes();
                    // Evaluate the basis functions at the evaluation nodes
                    const Eigen::MatrixXd basis =
                        sfl.EvalReferenceShapeFunctions(eval_nodes);
                    // Interpolate all basis functions
                    for (unsigned i = 0; i < basis.rows(); ++i) {
                      // Compute the basis function coefficients
                      const auto dofs = sfl.NodalValuesToDofs(basis.row(i));
                      // Only the i-th entry should be 1, all others 0
                      for (long j = 0; j < dofs.size(); ++j) {
                        ASSERT_NEAR(dofs[j], j == i ? 1 : 0, 1e-8)
                            << "o=" << to_char(o0) << to_char(o1) << to_char(o2)
                            << to_char(o3) << "p_interior=" << p_interior
                            << ", p_edges=" << p_edge0 << "," << p_edge1 << ","
                            << p_edge2 << "," << p_edge3 << " i=" << i
                            << " j=" << j << "\ndofs=[" << dofs << "]";
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, grad_segment) {
  // Test the gradients of the shape functions up to degree p=20
  quad::QuadRuleCache qr_cache;
  for (unsigned p = 1; p <= 20; ++p) {
    // Test for all possible combinations of orientations
    lf::mesh::Orientation orient[1];
    const lf::fe::FeHierarchicSegment<double> sfl(p, qr_cache);
    const auto eval_nodes = sfl.EvaluationNodes();
    // Compute the exact gradients
    const auto grad_exact = sfl.GradientsReferenceShapeFunctions(eval_nodes);
    // Compute the gradients with finite differences
    const double dx = std::sqrt(std::numeric_limits<double>::epsilon());
    const auto rsfp = sfl.EvalReferenceShapeFunctions(
        eval_nodes + Eigen::RowVectorXd::Constant(eval_nodes.cols(), dx / 2));
    const auto rsfn = sfl.EvalReferenceShapeFunctions(
        eval_nodes - Eigen::RowVectorXd::Constant(eval_nodes.cols(), dx / 2));
    const Eigen::MatrixXd grad_approx = (rsfp - rsfn) / dx;
    // Compare the gradients
    for (int i = 0; i < grad_exact.rows(); ++i) {
      const Eigen::RowVectorXd grad_i_exact = grad_exact.row(i);
      const Eigen::RowVectorXd grad_i_approx = grad_approx.row(i);
      const double max_diff =
          (grad_i_exact - grad_i_approx).array().abs().maxCoeff();
      ASSERT_TRUE(max_diff < 1e-8) << "i=" << i << std::endl;
    }
  }
}

TEST(fe_space_hierarchic, grad_tria) {
  // Test the gradients of the shape functions up to degree p=10
  quad::QuadRuleCache qr_cache;
  for (unsigned p_interior = 1; p_interior <= 6; ++p_interior) {
    for (unsigned p_edge0 = 1; p_edge0 <= 6; ++p_edge0) {
      for (unsigned p_edge1 = 1; p_edge1 <= 6; ++p_edge1) {
        for (unsigned p_edge2 = 1; p_edge2 <= 6; ++p_edge2) {
          // Test for all possible combinations of orientations
          lf::mesh::Orientation orient[4];
          for (const auto o0 : {lf::mesh::Orientation::positive,
                                lf::mesh::Orientation::negative}) {
            orient[0] = o0;
            for (const auto o1 : {lf::mesh::Orientation::positive,
                                  lf::mesh::Orientation::negative}) {
              orient[1] = o1;
              for (const auto o2 : {lf::mesh::Orientation::positive,
                                    lf::mesh::Orientation::negative}) {
                orient[2] = o2;
                for (const auto o3 : {lf::mesh::Orientation::positive,
                                      lf::mesh::Orientation::negative}) {
                  orient[3] = o3;
                  std::array<unsigned, 3> p_edges{{p_edge0, p_edge1, p_edge2}};
                  const lf::fe::FeHierarchicTria<double> sfl(
                      p_interior, p_edges, qr_cache, orient);
                  const auto eval_nodes = sfl.EvaluationNodes();
                  // Compute the exact gradients
                  const auto grad_exact =
                      sfl.GradientsReferenceShapeFunctions(eval_nodes);
                  // Compute the gradients with finite differences
                  const double d =
                      std::sqrt(std::numeric_limits<double>::epsilon());
                  Eigen::MatrixXd eval_nodes_xp = eval_nodes;
                  eval_nodes_xp.row(0) +=
                      Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                  Eigen::MatrixXd eval_nodes_xn = eval_nodes;
                  eval_nodes_xn.row(0) -=
                      Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                  Eigen::MatrixXd eval_nodes_yp = eval_nodes;
                  eval_nodes_yp.row(1) +=
                      Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                  Eigen::MatrixXd eval_nodes_yn = eval_nodes;
                  eval_nodes_yn.row(1) -=
                      Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                  const auto rsfxp =
                      sfl.EvalReferenceShapeFunctions(eval_nodes_xp);
                  const auto rsfxn =
                      sfl.EvalReferenceShapeFunctions(eval_nodes_xn);
                  const auto rsfyp =
                      sfl.EvalReferenceShapeFunctions(eval_nodes_yp);
                  const auto rsfyn =
                      sfl.EvalReferenceShapeFunctions(eval_nodes_yn);
                  const Eigen::MatrixXd dx_approx = (rsfxp - rsfxn) / d;
                  const Eigen::MatrixXd dy_approx = (rsfyp - rsfyn) / d;
                  // Compare the gradients
                  for (int i = 0; i < grad_exact.rows(); ++i) {
                    for (int j = 0; j < eval_nodes.cols(); ++j) {
                      const double dx_ij_exact = grad_exact(i, 2 * j + 0);
                      const double dy_ij_exact = grad_exact(i, 2 * j + 1);
                      const double dx_ij_approx = dx_approx(i, j);
                      const double dy_ij_approx = dy_approx(i, j);
                      ASSERT_TRUE(std::fabs(dx_ij_exact - dx_ij_approx) < 1e-5)
                          << "p=" << p_interior << " i=" << i << " eval_node=["
                          << eval_nodes.col(j).transpose() << "] orient=["
                          << static_cast<int>(orient[0]) << ", "
                          << static_cast<int>(orient[1]) << ", "
                          << static_cast<int>(orient[2]) << "]\ngrad_exact=["
                          << dx_ij_exact << ", " << dy_ij_exact
                          << "]\ngrad_approx=[" << dx_ij_approx << ", "
                          << dy_ij_approx << "]" << std::endl;
                      ASSERT_TRUE(std::fabs(dy_ij_exact - dy_ij_approx) < 1e-5)
                          << "p=" << p_interior << " i=" << i << " eval_node=["
                          << eval_nodes.col(j).transpose() << "] orient=["
                          << static_cast<int>(orient[0]) << ", "
                          << static_cast<int>(orient[1]) << ", "
                          << static_cast<int>(orient[2]) << "]\ngrad_exact=["
                          << dx_ij_exact << ", " << dy_ij_exact
                          << "]\ngrad_approx=[" << dx_ij_approx << ", "
                          << dy_ij_approx << "]" << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST(fe_space_hierarchic, grad_quad) {
  // Test the gradients of the shape functions up to degree p=6
  quad::QuadRuleCache qr_cache;
  unsigned max_p = 4;
#ifdef NDEBUG
  max_p = 6;
#endif

  for (unsigned p_interior = 1; p_interior <= max_p; ++p_interior) {
    for (unsigned p_edge0 = 1; p_edge0 <= max_p; ++p_edge0) {
      for (unsigned p_edge1 = 1; p_edge1 <= max_p; ++p_edge1) {
        for (unsigned p_edge2 = 1; p_edge2 <= max_p; ++p_edge2) {
          for (unsigned p_edge3 = 1; p_edge3 <= max_p; ++p_edge3) {
            // Test for all possible combinations of orientations
            lf::mesh::Orientation orient[4];
            for (const auto o0 : {lf::mesh::Orientation::positive,
                                  lf::mesh::Orientation::negative}) {
              orient[0] = o0;
              for (const auto o1 : {lf::mesh::Orientation::positive,
                                    lf::mesh::Orientation::negative}) {
                orient[1] = o1;
                for (const auto o2 : {lf::mesh::Orientation::positive,
                                      lf::mesh::Orientation::negative}) {
                  orient[2] = o2;
                  for (const auto o3 : {lf::mesh::Orientation::positive,
                                        lf::mesh::Orientation::negative}) {
                    orient[3] = o3;
                    std::array<unsigned, 4> edge_degrees{
                        {p_edge0, p_edge1, p_edge2, p_edge3}};
                    const lf::fe::FeHierarchicQuad<double> sfl(
                        p_interior, edge_degrees, qr_cache, orient);
                    const auto eval_nodes = sfl.EvaluationNodes();
                    // Compute the exact gradients
                    const auto grad_exact =
                        sfl.GradientsReferenceShapeFunctions(eval_nodes);
                    // Compute the gradients with finite differences
                    const double d =
                        std::sqrt(std::numeric_limits<double>::epsilon());
                    Eigen::MatrixXd eval_nodes_xp = eval_nodes;
                    eval_nodes_xp.row(0) +=
                        Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                    Eigen::MatrixXd eval_nodes_xn = eval_nodes;
                    eval_nodes_xn.row(0) -=
                        Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                    Eigen::MatrixXd eval_nodes_yp = eval_nodes;
                    eval_nodes_yp.row(1) +=
                        Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                    Eigen::MatrixXd eval_nodes_yn = eval_nodes;
                    eval_nodes_yn.row(1) -=
                        Eigen::RowVectorXd::Constant(eval_nodes.cols(), d / 2);
                    const auto rsfxp =
                        sfl.EvalReferenceShapeFunctions(eval_nodes_xp);
                    const auto rsfxn =
                        sfl.EvalReferenceShapeFunctions(eval_nodes_xn);
                    const auto rsfyp =
                        sfl.EvalReferenceShapeFunctions(eval_nodes_yp);
                    const auto rsfyn =
                        sfl.EvalReferenceShapeFunctions(eval_nodes_yn);
                    const Eigen::MatrixXd dx_approx = (rsfxp - rsfxn) / d;
                    const Eigen::MatrixXd dy_approx = (rsfyp - rsfyn) / d;
                    // Compare the gradients
                    for (int i = 0; i < grad_exact.rows(); ++i) {
                      for (int j = 0; j < eval_nodes.cols(); ++j) {
                        const double dx_ij_exact = grad_exact(i, 2 * j + 0);
                        const double dy_ij_exact = grad_exact(i, 2 * j + 1);
                        const double dx_ij_approx = dx_approx(i, j);
                        const double dy_ij_approx = dy_approx(i, j);
                        ASSERT_TRUE(std::fabs(dx_ij_exact - dx_ij_approx) <
                                    1e-8)
                            << "p="
                            << "p_interior=" << p_interior
                            << ", p_edges=" << p_edge0 << "," << p_edge1 << ","
                            << p_edge2 << "," << p_edge3 << " i=" << i
                            << " orient=[" << static_cast<int>(orient[0])
                            << ", " << static_cast<int>(orient[1]) << ", "
                            << static_cast<int>(orient[2]) << ", "
                            << static_cast<int>(orient[3]) << "]\ngrad_exact=["
                            << dx_ij_exact << ", " << dy_ij_exact
                            << "]\ngrad_approx=[" << dx_ij_approx << ", "
                            << dy_ij_approx << "]" << std::endl;
                        ASSERT_TRUE(std::fabs(dy_ij_exact - dy_ij_approx) <
                                    1e-8)
                            << "p="
                            << "p_interior=" << p_interior
                            << ", p_edges=" << p_edge0 << "," << p_edge1 << ","
                            << p_edge2 << "," << p_edge3 << " i=" << i
                            << " orient=[" << static_cast<int>(orient[0])
                            << ", " << static_cast<int>(orient[1]) << ", "
                            << static_cast<int>(orient[2]) << ", "
                            << static_cast<int>(orient[3]) << "]\ngrad_exact=["
                            << dx_ij_exact << ", " << dy_ij_exact
                            << "]\ngrad_approx=[" << dx_ij_approx << ", "
                            << dy_ij_approx << "]" << std::endl;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // end namespace lf::fe::test
