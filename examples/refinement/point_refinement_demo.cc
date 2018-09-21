/**
 * @file point_refinement_demo.cc
 * @brief example code for building and pointwise refining of a tensor product
 *        mesh
 *
 * This code generates a tensor-product mesh and performs a user specified
 * number of steps of pointwise uniform refinement on the cell containing the
 * point (0.5, 0.5).
 *
 * In the end the information about all meshes created in the process
 * is stored in the form of MATLAB functions.
 */

#include <iostream>

#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/refutils.h>
#include "lf/base/base.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"


int main() {
    using size_type = lf::base::size_type;
    using lf::mesh::utils::TikzOutputCtrl;

    std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
        std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);

    // Build single-cell tensor product mesh on unit square
    lf::mesh::hybrid2d::TPQuadMeshBuilder builder(mesh_factory_ptr);
    builder.setBottomLeftCorner(Eigen::Vector2d{0,0});
    builder.setTopRightCorner(Eigen::Vector2d{1,1});
    builder.setNoXCells(1);
    builder.setNoYCells(1);
    std::shared_ptr<lf::mesh::Mesh> mesh_ptr = builder.Build();

    // Output mesh information
    const lf::mesh::Mesh &mesh = *mesh_ptr;
    lf::mesh::utils::PrintInfo(mesh, std::cout);
    std::cout << std::endl;

    // Build mesh hierarchy
    lf::refinement::MeshHierarchy multi_mesh(mesh_ptr, mesh_factory_ptr);

    // Mark cell containing point (0.5, 0.5) which should be refined
    auto marker = [] (const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) {
        Eigen::MatrixXd ref_a(2, 1);
        Eigen::MatrixXd ref_b(2, 1);
        Eigen::MatrixXd ref_c(2, 1);
        Eigen::MatrixXd ref_d(2, 1);
        ref_a << 0,0;
        ref_b << 1,0;
        ref_c << 1,1;
        ref_d << 0,1;

        Eigen::MatrixXd point(2,1);
        point << .5,.5;

        for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
            auto geom = cell.Geometry();

            // Due to refinement we have to make sure that the global
            // coordinates are always ordered in the same way:
            // c ---- d
            //   |  |
            // a ---- b 
            std::vector<Eigen::MatrixXd> glob_coords = {geom->Global(ref_a), 
                geom->Global(ref_b), geom->Global(ref_c), geom->Global(ref_d)};

            auto coord_sort = [] (Eigen::MatrixXd x, Eigen::MatrixXd y) {
                return (x(1,0) < y(1,0)) || 
                       ((x(1,0) == y(1,0)) && (x(0,0) <= y(0,0)));
            };
            std::sort(glob_coords.begin(), glob_coords.end(), coord_sort);

            Eigen::MatrixXd glob_a = glob_coords[0];
            Eigen::MatrixXd glob_b = glob_coords[1];
            Eigen::MatrixXd glob_c = glob_coords[2];
            Eigen::MatrixXd glob_d = glob_coords[3];

            // Check if point lies inside cell
            if (glob_a(0,0) <= point(0,0) && glob_a(1,0) <= point(1,0) &&
                glob_b(0,0) >= point(0,0) && glob_b(1,0) <= point(1,0) &&
                glob_c(0,0) <= point(0,0) && glob_c(1,0) >= point(1,0) &&
                glob_d(0,0) >= point(0,0) && glob_d(1,0) >= point(1,0)) {

                // Check if cell contains edge
                auto edges = cell.SubEntities(1);
                if (edges.end() != std::find(edges.begin(),edges.end(),edge)) {
                    return true;
                }
            }
        }
        return false;
    };

    std::cout << "refinement steps: ";
    size_t refinement_steps;
    std::cin >> refinement_steps;

    std::cout << "pointwise refinement [0/1]: ";
    bool point_refinement;
    std::cin >> point_refinement;
    std::cout << std::endl;

    for (int step = 0; step < refinement_steps; ++step) {
        // Obtain pointer to mesh on finest level
        const size_type n_levels = multi_mesh.NumLevels();
        std::shared_ptr<const lf::mesh::Mesh> mesh_fine =
            multi_mesh.getMesh(n_levels - 1);

        // Print number of entities of various co-dimensions
        std::cout << "Mesh on level " << n_levels - 1 << ": " 
            << mesh_fine->Size(2)<< " nodes, " 
            << mesh_fine->Size(1) << " edges, " 
            << mesh_fine->Size(0) << " cells," << std::endl;

        lf::mesh::utils::writeTikZ(
                *mesh_fine, 
                std::string("refinement_mesh") + std::to_string(step) + ".txt",
                TikzOutputCtrl::RenderCells | TikzOutputCtrl::CellNumbering |
                TikzOutputCtrl::VerticeNumbering | 
                TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::EdgeNumbering);

        if (point_refinement) {
            multi_mesh.MarkEdges(marker);
            multi_mesh.RefineMarked();
        } else {
            multi_mesh.RefineRegular();
        }
    }

    // Generate  MATLAB functions that provide a description of all
    // levels of the mesh hierarchy
    std::cout << std::endl << "basename for MATLAB output: ";
    std::string basename;
    std::cin >> basename;
    lf::refinement::WriteMatlab(multi_mesh, basename);

    return 0;
}
