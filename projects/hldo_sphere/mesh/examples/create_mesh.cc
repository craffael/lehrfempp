#include <lf/io/vtk_writer.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <sphere_triag_mesh_builder.h>

#include <iostream>
#include <sstream>

/**
 * Creates example meshes with the refinement values
 *
 * creates meshes with refinement_level in smallest_mesh <= refinement_level <=
 * largest_mesh
 *
 * If verbose is the third input
 * print information about the smallest mesh specified
 * For each cell
 * local -> global mapping for edges
 * edge orientations for all edges
 *
 */
int main(int arg, char** args) {
  // prepare needed objects
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  // read out
  if (arg != 4 && arg != 5) {
    std::cerr << "Usage: " << args[0]
              << "smallest_mesh largest_mesh radius [verbose]" << '\n';
    std::quick_exit(1);
  }

  const double radius = std::stod(args[3]);
  sphere.setRadius(radius);

  const int largest_mesh = std::stoi(args[2]);
  const int smallest_mesh = std::stoi(args[1]);
  for (int i = 0; i <= largest_mesh; i++) {
    // create mesh
    // NOLINTNEXTLINE(clang-analyzer-deadcode.DeadStores)
    const int refinement_level = i + smallest_mesh;

    sphere.setRefinementLevel(i);
    const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();

    // check if verbositiy is activated
    // if so print mesh information of the smallest mesh
    if (i == 0 && arg == 5 && std::string(args[4]) == "verbose") {
      std::cout << "\n\n\nPoint information";
      for (const lf::mesh::Entity* point : mesh->Entities(2)) {
        const Eigen::MatrixXd vertex =
            lf::geometry::Corners(*(point->Geometry()));
        const lf::base::size_type index = mesh->Index(*point);
        std::cout << "\n\nPoint " << index << "\n" << vertex;
      }

      for (const lf::mesh::Entity* cell_p : mesh->Entities(0)) {
        const lf::base::size_type global_cell_index = mesh->Index(*cell_p);
        auto orientations = cell_p->RelativeOrientations();
        const Eigen::MatrixXd vertices =
            lf::geometry::Corners(*(cell_p->Geometry()));
        std::cout << "\n\n\nInfo for cell " << global_cell_index << "\n";
        std::cout << "Corners:\n" << vertices;

        int i = 0;
        for (const lf::mesh::Entity* edge_p : cell_p->SubEntities(1)) {
          const lf::base::size_type global_edge_index = mesh->Index(*edge_p);
          const int sign = lf::mesh::to_sign(orientations[i]);
          std::cout << "\n\nEdge local " << i << " global " << global_edge_index
                    << " orientation " << sign;
          i++;
        }
      }
    }

    // write mesh out with refinement_level and radius cut to 6 decimals
    const lf::io::VtkWriter vtk_writer_1(
        mesh, "meshs/sphere_" + std::to_string(i) + "_" +
                  std::to_string(radius) + ".vtk");
  }
}
