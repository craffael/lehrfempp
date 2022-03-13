#define _USE_MATH_DEFINES
#include <lf/io/vtk_writer.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <sphere_triag_mesh_builder.h>

/**
 * Creates example meshes with the refinement values
 * 0
 * 1
 * 2
 * 3
 * 4
 * 5
 * 6
 * 7
 */
int main(int arg, char** args) {
  // prepare needed objects
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  sphere.setRadius(1);

  sphere.setRefinementLevel(0);
  const std::shared_ptr<lf::mesh::Mesh> mesh_0 = sphere.Build();
  lf::io::VtkWriter vtk_writer_0(mesh_0, "sphere_0.vtk");

  sphere.setRefinementLevel(1);
  const std::shared_ptr<lf::mesh::Mesh> mesh_1 = sphere.Build();
  lf::io::VtkWriter vtk_writer_1(mesh_1, "sphere_1.vtk");

  sphere.setRefinementLevel(2);
  const std::shared_ptr<lf::mesh::Mesh> mesh_2 = sphere.Build();
  lf::io::VtkWriter vtk_writer_2(mesh_2, "sphere_2.vtk");

  sphere.setRefinementLevel(3);
  const std::shared_ptr<lf::mesh::Mesh> mesh_3 = sphere.Build();
  lf::io::VtkWriter vtk_writer_3(mesh_3, "sphere_3.vtk");

  sphere.setRefinementLevel(4);
  const std::shared_ptr<lf::mesh::Mesh> mesh_4 = sphere.Build();
  lf::io::VtkWriter vtk_writer_4(mesh_4, "sphere_4.vtk");

  sphere.setRefinementLevel(5);
  const std::shared_ptr<lf::mesh::Mesh> mesh_5 = sphere.Build();
  lf::io::VtkWriter vtk_writer_5(mesh_5, "sphere_5.vtk");

  sphere.setRefinementLevel(6);
  const std::shared_ptr<lf::mesh::Mesh> mesh_6 = sphere.Build();
  lf::io::VtkWriter vtk_writer_6(mesh_6, "sphere_6.vtk");

  sphere.setRefinementLevel(7);
  const std::shared_ptr<lf::mesh::Mesh> mesh_7 = sphere.Build();
  lf::io::VtkWriter vtk_writer_7(mesh_7, "sphere_7.vtk");
}
