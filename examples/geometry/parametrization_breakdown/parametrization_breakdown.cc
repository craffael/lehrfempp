/**
 * @file
 * @brief Investigates when parametrization from reference element breaks down
 * @author Anian Ruoss
 * @date   2019-02-04 13:36:17
 * @copyright MIT License
 */

#include <lf/geometry/geometry.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

#include <Eigen/Eigen>
#include <filesystem>
#include <fstream>
#include <string>

/**
 * @brief Stores an Eigen::MatrixXd to .csv file
 * @param file_path path to .csv file
 * @param matrix matrix to be saved
 */
void writeMatrixToCSV(const std::string& file_path,
                      const Eigen::MatrixXd& matrix) {
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file(file_path.c_str());
  file << matrix.format(CSVFormat);
  file.close();
}

/**
 * @brief Saves the global vertex/midpoint coordinates
 * @param file_path path to .csv file
 * @param geom second-order geometry object
 */
void storeSecondOrderCoords(const std::string& file_path,
                            const lf::geometry::Geometry& geom) {
  const Eigen::MatrixXd& local_vertex_coords = geom.RefEl().NodeCoords();
  const long num_vertices = local_vertex_coords.cols();

  Eigen::MatrixXd all_local_coords(local_vertex_coords.rows(),
                                   2 * num_vertices);

  // compute local vertex and midpoint coordinates from reference element
  for (auto col_idx = 0; col_idx < num_vertices; ++col_idx) {
    all_local_coords.col(col_idx) = local_vertex_coords.col(col_idx);
    all_local_coords.col(col_idx + num_vertices) =
        (local_vertex_coords.col(col_idx) +
         local_vertex_coords.col((col_idx + 1) % num_vertices)) /
        2.;
  }

  writeMatrixToCSV(file_path, geom.Global(all_local_coords));
}

/**
 * @brief Evaluates parametrization and corresponding jacobian determinant
 * @param base_file_path base path for .csv files
 * @param geom geometry object
 * @param qr_order order of quadrature rule to sample points from
 */
void storeParametrizationEvals(const std::string& base_file_path,
                               const lf::geometry::Geometry& geom,
                               const size_t& qr_order) {
  // use quadrature to sample random points on the reference element
  auto qr = lf::quad::make_QuadRule(geom.RefEl(), qr_order);

  const auto& points = qr.Points();
  const Eigen::MatrixXd& jacobians = geom.Jacobian(points);
  Eigen::VectorXd determinants(points.cols());

  for (Eigen::Index point_idx = 0; point_idx < points.cols(); ++point_idx) {
    determinants(point_idx) = jacobians
                                  .block(0, point_idx * geom.DimLocal(),
                                         geom.DimGlobal(), geom.DimLocal())
                                  .determinant();
  }

  writeMatrixToCSV(base_file_path + "_refpoints.csv", points);
  writeMatrixToCSV(base_file_path + "_points.csv", geom.Global(points));
  writeMatrixToCSV(base_file_path + "_jacdets.csv", determinants);
}

/**
 * @brief Computes volume of geometry object by means of overkill quadrature
 * @param geom geometry object
 * @return geometry volume
 */
double computeGeometryVolume(const lf::geometry::Geometry& geom) {
  const auto qr = lf::quad::make_QuadRule(geom.RefEl(), 23);

  const auto& points = qr.Points();
  const auto& weights = qr.Weights();
  const auto& integrationElements = geom.IntegrationElement(points);

  double vol = 0.;

  for (Eigen::Index j = 0; j < points.cols(); ++j) {
    vol += weights(j) * integrationElements(j);
  }

  return vol;
}

int main() {
  // create a directory to store results
  const std::string results_dir = "results/";
  std::filesystem::create_directory(results_dir);

  // define second-order geometry elements
  lf::geometry::TriaO2 tria(
      (Eigen::MatrixXd(2, 6) << 1, 6, 3, 3.7, 4.2, 2.3, 1, 3, 8, 1.2, 5.2, 4.5)
          .finished());
  lf::geometry::TriaO2 tria_degenerate(
      (Eigen::MatrixXd(2, 6) << 1, 6, 3, 5, 4.5, 1.75, 1, 3, 8, 4, 9, 6.5)
          .finished());
  lf::geometry::QuadO2 quad((Eigen::MatrixXd(2, 8) << 3, 7, 4, 1, 5, 5.4, 2.5,
                             1.8, 1, 3, 7, 8, 2.5, 5.9, 6.9, 4.1)
                                .finished());
  lf::geometry::QuadO2 quad_degenerate(
      (Eigen::MatrixXd(2, 8) << 3, 7, 4, 1, 6, 5, 2, 2, 1, 3, 7, 8, 0, 6, 5, 2)
          .finished());

  // store vertex/midpoint coordinates, random point evaluations and
  // corresponding jacobian determinants for every geometry object
  for (const auto& geom_element :
       {std::pair<std::string, lf::geometry::Geometry*>{"tria", &tria},
        std::pair<std::string, lf::geometry::Geometry*>{"tria_degenerate",
                                                        &tria_degenerate},
        std::pair<std::string, lf::geometry::Geometry*>{"quad", &quad},
        std::pair<std::string, lf::geometry::Geometry*>{"quad_degenerate",
                                                        &quad_degenerate}}) {
    storeSecondOrderCoords(results_dir + geom_element.first + "_coords.csv",
                           *geom_element.second);
    storeParametrizationEvals(results_dir + geom_element.first,
                              *geom_element.second, 50);

    const double volume = computeGeometryVolume(*geom_element.second);
    double refined_volume = 0;

    // compute child geometries by means of regular refinement
    auto children = geom_element.second->ChildGeometry(
        lf::refinement::Hybrid2DRefinementPattern(geom_element.second->RefEl(),
                                                  lf::refinement::rp_regular),
        0);

    // store vertex/midpoint coordinates, random point evaluations and
    // corresponding jacobian determinants for every child geometry object
    for (std::size_t child_idx = 0; child_idx < children.size(); ++child_idx) {
      storeSecondOrderCoords(results_dir + geom_element.first + "_child_" +
                                 std::to_string(child_idx) + "_coords.csv",
                             *children[child_idx]);
      storeParametrizationEvals(results_dir + geom_element.first + "_child_" +
                                    std::to_string(child_idx),
                                *children[child_idx], 20);

      refined_volume += computeGeometryVolume(*children[child_idx]);
    }

    // save volumes
    writeMatrixToCSV(results_dir + geom_element.first + "_volumes.csv",
                     (Eigen::VectorXd(2) << volume, refined_volume).finished());
  }

  return 0;
}
