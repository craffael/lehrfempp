#include "dirac_convergence_test.h"

#include "debugging.h"

namespace projects::hldo_sphere::debugging {

/**
 *
 * @brief Solves the hodge laplacian source problems for the tensor product of
 * passed refinement levels and ks
 *
 * @param refinement_levels integer list containig all the levels
 * @param ks list of all ks to be used
 *
 */
void DiracConvergenceTest::Compute(unsigned refinement_levels) {
  // Initialize solution wrapper
  SolutionList solutions;
  std::vector<std::shared_ptr<const lf::mesh::Mesh>> meshs =
      std::vector<std::shared_ptr<const lf::mesh::Mesh>>(refinement_levels);

  projects::hldo_sphere::operators::DiracOperatorSourceProblem lse_builder;

  // functions are passed by reference hence changing the k still influences
  // the functions
  lse_builder.SetLoadFunctions(f_zero_, f_one_, f_two_);
  lse_builder.SetK(k_);

  // Read the mesh from the gmsh file
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  // we always take the same radius
  sphere.setRadius(1.);

  // start timer for total time
  const std::chrono::steady_clock::time_point start_time_total =
      std::chrono::steady_clock::now();

  // Initialize refinement level
  solutions.mu_zero.resize(refinement_levels);
  solutions.mu_one.resize(refinement_levels);
  solutions.mu_two.resize(refinement_levels);

  for (unsigned il = 0; il < refinement_levels; ++il) {
    std::cout << "\nStart computation of refinement_level " << il << std::flush;

    // start timer
    std::chrono::steady_clock::time_point start_time =
        std::chrono::steady_clock::now();

    // get mesh
    sphere.setRefinementLevel(il);
    const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();
    meshs[il] = mesh;

    // end timer
    std::chrono::steady_clock::time_point end_time =
        std::chrono::steady_clock::now();
    auto elapsed_sec =
        lf::base::narrow<double>(
            std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                  start_time)
                .count()) /
        1000.;

    std::cout << " -> Built Mesh " << elapsed_sec << " [s] " << std::flush;

    // setup the problem
    lse_builder.SetMesh(mesh);

    // start timer
    start_time = std::chrono::steady_clock::now();

    // compute the system
    lse_builder.Compute();

    // end timer
    end_time = std::chrono::steady_clock::now();
    elapsed_sec = lf::base::narrow<double>(
                      std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_time - start_time)
                          .count()) /
                  1000.;

    std::cout << " -> Computed LSE " << elapsed_sec << " [s] " << std::flush;

    // start timer
    start_time = std::chrono::steady_clock::now();

    // solve the system
    lse_builder.Solve();

    // end timer
    end_time = std::chrono::steady_clock::now();

    elapsed_sec = lf::base::narrow<double>(
                      std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_time - start_time)
                          .count()) /
                  1000.;

    std::cout << " -> Solved System " << elapsed_sec << " [s] " << std::flush;

    // store solutions
    solutions.mu_zero[il] = lse_builder.GetMu(0);
    solutions.mu_one[il] = lse_builder.GetMu(1);
    solutions.mu_two[il] = lse_builder.GetMu(2);
  }  // end loop level

  // end timer total
  const std::chrono::steady_clock::time_point end_time_total =
      std::chrono::steady_clock::now();
  const double elapsed_sec_total =
      lf::base::narrow<double>(
          std::chrono::duration_cast<std::chrono::milliseconds>(
              end_time_total - start_time_total)
              .count()) /
      1000.;

  std::cout << "\nTotal computation time for all levels " << elapsed_sec_total
            << " [s]\n";

  /* Do convergence analyisis of the solutions that is
   * | \|u_k\| - \|u_{k-1}\| |^2 -> 0
   * with algebraic rate is expected
   */

  // create results direcotry
  const std::string main_dir = "results";
  std::filesystem::create_directory(main_dir);

  // prepare output
  const int table_width = 13;
  const int precision = 4;

  // create csv file
  std::ofstream csv_file;
  std::string csv_name = concat("result_dirac_convergence_", refinement_levels);
  std::replace(csv_name.begin(), csv_name.end(), '.', '_');
  csv_file.open(concat(main_dir, "/", csv_name, ".csv"));
  csv_file << "numCells,"
           << "numEdges,"
           << "numVerts,"
           << "hMax";

  // define square functions for norms
  auto square_scalar = [](complex a) -> double {
    return std::abs(a) * std::abs(a);
  };
  auto square_vector =
      [](const Eigen::Matrix<complex, Eigen::Dynamic, 1> &a) -> double {
    return a.squaredNorm();
  };

  // define quadrule for norms
  const lf::quad::QuadRule qr =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 4);

  Eigen::VectorXd L2ErrorZero = Eigen::VectorXd::Zero(refinement_levels);
  Eigen::VectorXd L2ErrorOne = Eigen::VectorXd::Zero(refinement_levels);
  Eigen::VectorXd L2ErrorTwo = Eigen::VectorXd::Zero(refinement_levels);
  Eigen::VectorXd hMax = Eigen::VectorXd::Zero(refinement_levels);

  // introduce columns for errors
  csv_file << ",L2ErrorZero,DiffOfL2Zero,L2ErrorOne,DiffOfL2One,"
              "L2ErrorTwo,DiffOfL2Two";
  csv_file << '\n';

  // loop over all levels contained in the solution
  for (lf::base::size_type lvl = 0; lvl <= refinement_levels; ++lvl) {
    // create mesh Functions for solutions the level
    auto &sol_mesh = meshs[lvl];

    // compute meshwidth
    for (const lf::mesh::Entity *e : sol_mesh->Entities(1)) {
      const double h = lf::geometry::Volume(*(e->Geometry()));
      if (h > hMax(lvl)) {
        hMax(lvl) = h;
      }
    }

    // print level and mesh informations
    csv_file << sol_mesh->NumEntities(0) << "," << sol_mesh->NumEntities(1)
             << "," << sol_mesh->NumEntities(2) << "," << hMax(lvl);

    const projects::hldo_sphere::post_processing::MeshFunctionWhitneyZero
        mf_zero(solutions.mu_zero[lvl], sol_mesh);
    const projects::hldo_sphere::post_processing::MeshFunctionWhitneyOne mf_one(
        solutions.mu_one[lvl], sol_mesh);
    const projects::hldo_sphere::post_processing::MeshFunctionWhitneyTwo mf_two(
        solutions.mu_two[lvl], sol_mesh);

    // get error for zero form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_zero =
        projects::hldo_sphere::post_processing::L2norm(sol_mesh, mf_zero,
                                                       square_scalar, qr);
    L2ErrorZero(lvl) = std::get<0>(L2_zero);

    // get error for one form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_one =
        projects::hldo_sphere::post_processing::L2norm(sol_mesh, mf_one,
                                                       square_vector, qr);
    L2ErrorOne(lvl) = std::get<0>(L2_one);

    // get error for two form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_two =
        projects::hldo_sphere::post_processing::L2norm(sol_mesh, mf_two,
                                                       square_scalar, qr);
    L2ErrorTwo(lvl) = std::get<0>(L2_two);

    double l2RateZero = 0;
    double l2RateOne = 0;
    double l2RateTwo = 0;

    // we can only compute order form the second level
    if (lvl > 0) {
      l2RateZero = abs(L2ErrorZero(lvl) - L2ErrorZero(lvl - 1));
      l2RateOne = abs(L2ErrorOne(lvl) - L2ErrorOne(lvl - 1));
      l2RateTwo = abs(L2ErrorTwo(lvl) - L2ErrorTwo(lvl - 1));
    }

    /******************************
     * append solution of the current k in the outputs
     ******************************/

    csv_file << "," << L2ErrorZero(lvl) << "," << l2RateZero << ","
             << L2ErrorOne(lvl) << "," << l2RateOne << "," << L2ErrorTwo(lvl)
             << "," << l2RateTwo;
    csv_file << "\n";
  }  // end loop over levels

  // close csv file
  csv_file.close();
};

}  // namespace projects::hldo_sphere::debugging
