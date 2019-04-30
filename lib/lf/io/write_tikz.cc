/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/base/base.h"
#include "lf/mesh/mesh.h"

#include <fstream>

namespace lf::io {

TikzOutputCtrl operator|(const TikzOutputCtrl &lhs, const TikzOutputCtrl &rhs) {
  return static_cast<TikzOutputCtrl>(static_cast<unsigned int>(lhs) |
                                     static_cast<unsigned int>(rhs));
}

TikzOutputCtrl operator&(const TikzOutputCtrl &lhs, const TikzOutputCtrl &rhs) {
  return static_cast<TikzOutputCtrl>(static_cast<unsigned int>(lhs) &
                                     static_cast<unsigned int>(rhs));
}

bool ControlPointsCubicBezier(const Eigen::Matrix2d &vertices,
                              const Eigen::Vector2d &midpoint,
                              Eigen::Matrix2d &control_points) {
  auto collinear = [](const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                      const Eigen::Vector2d &p2) {
    Eigen::Matrix2d tmp;
    tmp << p1 - p0, p2 - p0;

    return std::abs(tmp.determinant()) < 1e-9;
  };

  if (collinear(vertices.col(0), midpoint, vertices.col(1))) {
    return true;
  }

  const double x1 = vertices(0, 0);
  const double y1 = vertices(1, 0);
  const double x2 = vertices(0, 1);
  const double y2 = vertices(1, 1);
  const double x3 = midpoint(0);
  const double y3 = midpoint(1);

  // compute coefficients of 2D parabola through points
  const double denominator = (x1 - x2) * (x1 - x3) * (x2 - x3);

  const double alpha =
      (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denominator;
  const double beta =
      (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) /
      denominator;

  // define derivative of parabola
  auto df = [&alpha, &beta](double x) { return 2 * alpha * x + beta; };

  // compute control point for the quadratic Bézier curve
  Eigen::Vector2d quadratic_control_point;
  quadratic_control_point << (x1 + x2) / 2., y1 + df(x1) * (x2 - x1) / 2.;

  // compute control points for cubic Bézier curve
  control_points << (2. * quadratic_control_point + vertices.col(0)) / 3.,
      (2. * quadratic_control_point + vertices.col(1)) / 3.;

  return false;
}

bool writeTikZ(const lf::mesh::Mesh &mesh, const std::string &filename,
               std::function<bool(const lf::mesh::Entity &)> &&selector,
               TikzOutputCtrl output_ctrl) {
  // Open output file for writing and quit in case of failure
  std::ofstream outfile(filename);
  if (!outfile) {
    return false;
  }
  // ----------------------------------------------------------------
  // For the enum flags: TikzOutputCtrl
  bool EdgeNumOn =
      static_cast<bool>(output_ctrl & TikzOutputCtrl::EdgeNumbering);
  bool NodeNumOn =
      static_cast<bool>(output_ctrl & TikzOutputCtrl::NodeNumbering);
  bool CellNumOn =
      static_cast<bool>(output_ctrl & TikzOutputCtrl::CellNumbering);
  bool VerticeNumOn =
      static_cast<bool>(output_ctrl & TikzOutputCtrl::VerticeNumbering);
  bool RenderCellsOn =
      static_cast<bool>(output_ctrl & TikzOutputCtrl::RenderCells);
  bool ArrowsOn = static_cast<bool>(output_ctrl & TikzOutputCtrl::ArrowTips);
  bool SecondOrder =
      static_cast<bool>(output_ctrl & TikzOutputCtrl::SecondOrder);

  using size_type = std::size_t;         // lf::base::size_type;
  using dim_t = lf::base::RefEl::dim_t;  // lf::base::dim_t;
  const Eigen::MatrixXd zero(Eigen::MatrixXd::Zero(0, 1));

  // Obtain topological dimension of the mesh
  const dim_t dim_mesh = mesh.DimMesh();
  LF_VERIFY_MSG(dim_mesh == 2, "writeTikZ() only available for 2D meshes");

  // Run through nodes
  const dim_t node_codim(dim_mesh);  // Codimension number for nodes in the mesh
  const size_type no_of_nodes =
      mesh.NumEntities(node_codim);  // No. of nodes in codimension for nodes
  size_type node_count = 0;

  // START writing to file
  outfile << "% TikZ document graphics \n";

  // Scale font size for large meshes
  if (no_of_nodes >= 25) {
    outfile << "\\begin{tikzpicture}[scale=6, >= stealth, inner sep=0pt, "
               "minimum size=0.2cm]\n";
    outfile << "\\tikzstyle{every node}=[font=\\tiny]\n";
  } else {
    outfile << "\\begin{tikzpicture}[scale=4, >= stealth, inner sep=0pt, "
               "minimum size=0.35cm]\n";
  }

  // ---------------------------------------------------------------------------

  // Loop through codimensions
  for (int co_dim = 0; co_dim <= dim_mesh; co_dim++) {
    // Loop through all types of entities
    for (const mesh::Entity &obj : mesh.Entities(co_dim)) {
      if (selector(obj)) {  // IF SELECTOR -------
        size_type obj_idx = mesh.Index(obj);
        lf::base::RefEl obj_refel = obj.RefEl();
        int num_nodes_obj = obj_refel.NumNodes();
        const geometry::Geometry *obj_geo_ptr = obj.Geometry();
        const Eigen::MatrixXd &obj_corners(obj_refel.NodeCoords());
        const Eigen::MatrixXd vertices = obj_geo_ptr->Global(obj_corners);
        const Eigen::MatrixXd center =
            vertices.rowwise().sum() / vertices.cols();
        Eigen::MatrixXd center_mat(center.rows(), num_nodes_obj);

        switch (obj_refel) {
          case lf::base::RefEl::kPoint(): {
            if (NodeNumOn) {
              outfile << "\\draw[red, fill = white] (" << vertices(0, 0) << ","
                      << vertices(1, 0) << ") "
                      << "node[draw, circle, fill = white] {" << obj_idx
                      << "};\n";
            } else {
              outfile << "\\draw[red] (" << vertices(0, 0) << ","
                      << vertices(1, 0) << ") "
                      << "node[] {*};\n";
            }  // if NodeNumOn
            break;

          }  // case kPoint

          case lf::base::RefEl::kSegment(): {
            center_mat << center, center;
            const Eigen::MatrixXd scaled_vertices =
                vertices * 0.80 + center_mat * 0.2;
            const Eigen::MatrixXd semi_scaled_vertices =
                vertices * 0.95 + center_mat * 0.05;

            if (ArrowsOn) {
              if (NodeNumOn) {
                outfile << "\\draw[->] (" << scaled_vertices(0, 0) << ","
                        << scaled_vertices(1, 0) << ")";
              } else {
                outfile << "\\draw[->] (" << semi_scaled_vertices(0, 0) << ","
                        << semi_scaled_vertices(1, 0) << ")";
              }

              std::string edge = " -- ";
              std::string midpoint_node =
                  " node[pos=0.5,black] {" + std::to_string(obj_idx) + "}";

              if (SecondOrder) {
                Eigen::Vector2d midpoint = obj_geo_ptr->Global(
                    (Eigen::MatrixXd(1, 1) << 0.5).finished());

                Eigen::Matrix2d control;
                bool collinear =
                    ControlPointsCubicBezier(vertices, midpoint, control);

                if (!collinear) {
                  edge = " .. controls (" + std::to_string(control(0, 0)) +
                         "," + std::to_string(control(1, 0)) + ")" + " and " +
                         "(" + std::to_string(control(0, 1)) + "," +
                         std::to_string(control(1, 1)) + ")" + " .. ";
                  midpoint_node =
                      "; \\draw (" + std::to_string(midpoint(0, 0)) + "," +
                      std::to_string(midpoint(1, 0)) + ")" + " node[black] {" +
                      std::to_string(obj_idx) + "}";
                }
              }

              outfile << edge;

              if (NodeNumOn) {
                outfile << "(" << scaled_vertices(0, 1) << ","
                        << scaled_vertices(1, 1) << ")";
              } else {
                outfile << "(" << semi_scaled_vertices(0, 1) << ","
                        << semi_scaled_vertices(1, 1) << ")";
              }

              if (EdgeNumOn) {
                outfile << midpoint_node;
              }

              outfile << ";\n";

            } else {
              outfile << "\\draw[] "
                      << "(" << vertices(0, 0) << "," << vertices(1, 0) << ")";

              std::string edge = " -- ";
              std::string midpoint_node =
                  " node[pos=0.5,black] {" + std::to_string(obj_idx) + "}";

              if (SecondOrder) {
                Eigen::Vector2d midpoint = obj_geo_ptr->Global(
                    (Eigen::MatrixXd(1, 1) << 0.5).finished());

                Eigen::Matrix2d control;
                bool collinear =
                    ControlPointsCubicBezier(vertices, midpoint, control);

                if (!collinear) {
                  edge = " .. controls (" + std::to_string(control(0, 0)) +
                         "," + std::to_string(control(1, 0)) + ")" + " and " +
                         "(" + std::to_string(control(0, 1)) + "," +
                         std::to_string(control(1, 1)) + ")" + " .. ";
                  midpoint_node = " (" + std::to_string(midpoint(0, 0)) + "," +
                                  std::to_string(midpoint(1, 0)) + ")" +
                                  " node[black] {" + std::to_string(obj_idx) +
                                  "}";
                }
              }

              outfile << edge;

              outfile << "(" << vertices(0, 1) << "," << vertices(1, 1) << ")";

              if (EdgeNumOn) {
                outfile << midpoint_node;
              }

              outfile << ";\n";
            }  // arrows on

            break;
          }  // case kSegment

          case lf::base::RefEl::kTria(): {
            center_mat << center, center, center;
            const Eigen::MatrixXd scaled_vertices =
                vertices * 0.70 + center_mat * 0.30;

            if (RenderCellsOn) {
              if (VerticeNumOn) {
                outfile << "\\draw[green] (" << scaled_vertices(0, 0) << ","
                        << scaled_vertices(1, 0) << ") node[] {0} -- ("
                        << scaled_vertices(0, 1) << "," << scaled_vertices(1, 1)
                        << ") node[] {1} -- (" << scaled_vertices(0, 2) << ","
                        << scaled_vertices(1, 2) << ") node[] {2} -- cycle;\n";
              } else {
                outfile << "\\draw[green] (" << scaled_vertices(0, 0) << ","
                        << scaled_vertices(1, 0) << ") -- ("
                        << scaled_vertices(0, 1) << "," << scaled_vertices(1, 1)
                        << ") -- (" << scaled_vertices(0, 2) << ","
                        << scaled_vertices(1, 2) << ") -- cycle;\n";
              }  // if EdgeNumOn

              if (CellNumOn) {
                outfile << "\\draw[green] (" << center(0, 0) << ","
                        << center(1, 0) << ") node[] {" << obj_idx << "};\n";
              }

            }  // RenderCellsOn
            break;
          }  // case kTria

          case lf::base::RefEl::kQuad(): {
            center_mat << center, center, center, center;
            const Eigen::MatrixXd scaled_vertices =
                vertices * 0.70 + center_mat * 0.3;

            if (RenderCellsOn) {
              if (VerticeNumOn) {
                outfile << "\\draw[magenta] (" << scaled_vertices(0, 0) << ","
                        << scaled_vertices(1, 0) << ") node[] {0} -- ("
                        << scaled_vertices(0, 1) << "," << scaled_vertices(1, 1)
                        << ") node[] {1} -- (" << scaled_vertices(0, 2) << ","
                        << scaled_vertices(1, 2) << ") node[] {2} -- ("
                        << scaled_vertices(0, 3) << "," << scaled_vertices(1, 3)
                        << ") node[] {3} -- cycle;\n";
              } else {
                outfile << "\\draw[magenta] (" << scaled_vertices(0, 0) << ","
                        << scaled_vertices(1, 0) << ") -- ("
                        << scaled_vertices(0, 1) << "," << scaled_vertices(1, 1)
                        << ") -- (" << scaled_vertices(0, 2) << ","
                        << scaled_vertices(1, 2) << ") -- ("
                        << scaled_vertices(0, 3) << "," << scaled_vertices(1, 3)
                        << ") -- cycle;\n";
              }

              if (CellNumOn) {
                outfile << "\\draw[magenta] (" << center(0, 0) << ","
                        << center(1, 0) << ") node[] {" << obj_idx << "};\n";
              }

            }  // RenderCellsOn

            break;
          }  // case kQuad

          default: {
            std::cout << "Error for object " << obj_idx << " in co-dim "
                      << co_dim << std::endl;
            std::cout << "Object type: " << obj_refel << std::endl;
            break;
          }  // default
        }    // switch
      }      // IF SELECTOR --------------------
    }        // for entities

    node_count++;
    // LF_VERIFY_MSG(node_count == no_of_nodes, "Node count mismatch");
  }  // for codim

  outfile << "\\end{tikzpicture}\n";
  return true;
}  // writetikz

// Second version of writeTikZ using default selector
bool writeTikZ(const lf::mesh::Mesh &mesh, const std::string &filename,
               TikzOutputCtrl output_ctrl) {
  return writeTikZ(mesh, filename,
                   [](const lf::mesh::Entity &) { return true; }, output_ctrl);
}

}  // namespace lf::io
