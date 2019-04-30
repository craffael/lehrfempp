/**
 * @file
 * @brief Implementation of writeMatplotlib
 * @author Anian Ruoss
 * @date   2018-10-08 18:27:17
 * @copyright MIT License
 */

#include "write_matplotlib.h"

#include <fstream>

namespace lf::io {

void writeMatplotlib(const lf::mesh::Mesh &mesh, std::string filename,
                     bool second_order) {
  using dim_t = lf::base::RefEl::dim_t;

  // append .txt to filename if necessary
  if (filename.compare(filename.size() - 4, 4, ".txt") != 0) {
    filename += ".txt";
  }

  std::ofstream file(filename);

  if (file.is_open()) {
    const dim_t dim_mesh = mesh.DimMesh();
    LF_VERIFY_MSG(dim_mesh == 2,
                  "write_matplotlib() only available for 2D meshes");

    // store points to file
    {
      std::vector<std::pair<size_t, Eigen::Vector2d>> points;

      for (const lf::mesh::Entity &point : mesh.Entities(2)) {
        points.emplace_back(std::make_pair(
            mesh.Index(point),
            point.Geometry()->Global(point.RefEl().NodeCoords())));
      }

      for (const auto &point : points) {
        const Eigen::Vector2d coords = point.second;
        file << "Point"
             << " " << point.first << " " << coords(0) << " " << coords(1)
             << std::endl;
      }
    }

    // store segments to file
    {
      std::vector<std::pair<size_t, Eigen::Matrix<double, 2, Eigen::Dynamic>>>
          segments;

      for (const lf::mesh::Entity &segment : mesh.Entities(1)) {
        if (second_order) {
          segments.emplace_back(std::make_pair(
              mesh.Index(segment),
              segment.Geometry()->Global(
                  (Eigen::MatrixXd(1, 3) << 0, 1, 0.5).finished())));
        } else {
          segments.emplace_back(
              std::make_pair(mesh.Index(segment),
                             segment.Geometry()->Global(
                                 (Eigen::MatrixXd(1, 2) << 0, 1).finished())));
        }
      }

      for (const auto &segment : segments) {
        file << (second_order ? "SegmentO2" : "SegmentO1") << " "
             << segment.first;
        const Eigen::Matrix<double, 2, Eigen::Dynamic> coords = segment.second;

        for (int col = 0; col < coords.cols(); ++col) {
          file << " " << coords(0, col) << " " << coords(1, col);
        }

        file << std::endl;
      }
    }

    {
      std::vector<std::pair<size_t, Eigen::Matrix<double, 2, Eigen::Dynamic>>>
          cells;

      for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
        const lf::geometry::Geometry *geometry_ptr = cell.Geometry();
        const Eigen::MatrixXd local_vertices = cell.RefEl().NodeCoords();
        const Eigen::MatrixXd vertices = geometry_ptr->Global(local_vertices);

        if (!second_order) {
          cells.emplace_back(std::make_pair(mesh.Index(cell), vertices));

        } else {
          // compute midpoints
          const int num_points = vertices.cols();
          Eigen::MatrixXd lower_diagonal =
              Eigen::MatrixXd::Zero(num_points, num_points);
          lower_diagonal.diagonal<-1>() = Eigen::VectorXd::Ones(num_points - 1);
          lower_diagonal(0, num_points - 1) = 1;

          const Eigen::MatrixXd local_midpoints =
              local_vertices *
              (.5 * (lower_diagonal +
                     Eigen::MatrixXd::Identity(num_points, num_points)));
          const Eigen::MatrixXd midpoints =
              geometry_ptr->Global(local_midpoints);

          cells.emplace_back(std::make_pair(
              mesh.Index(cell),
              (Eigen::MatrixXd(2, 2 * num_points) << vertices, midpoints)
                  .finished()));
        }
      }

      for (const auto &cell : cells) {
        const Eigen::Matrix<double, 2, Eigen::Dynamic> coords = cell.second;

        if (coords.cols() == 3 || coords.cols() == 6) {
          file << (second_order ? "TriaO2" : "TriaO1");
        } else if (coords.cols() == 4 || coords.cols() == 8) {
          file << (second_order ? "QuadO2" : "QuadO1");

        } else {
          LF_VERIFY_MSG(false, "Unknown cell geometry");
        }

        file << " " << cell.first;

        for (int col = 0; col < coords.cols(); ++col) {
          file << " " << coords(0, col) << " " << coords(1, col);
        }

        file << std::endl;
      }
    }
  }
}

}  // namespace lf::io
