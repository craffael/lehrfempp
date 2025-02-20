// This file contains Doxygen snippets for the quick reference mesh document
// It defines a function to hold all code snippets and includes necessary
// imports.

#include <lf/io/io.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

#include <iostream>
#include <memory>
#include <span>

void qr_mesh_snippets() {
  {
    //! [mesh Overview]
    // Generate a simple test mesh
    const std::shared_ptr<const lf::mesh::Mesh> mesh_ptr =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

    // Auto can be used to simplify the type declaration
    auto mesh_ptr_auto = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
    //! [mesh Overview]
  }

  //! [mesh Mesh Access]
  // create a pointer to the mesh
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
  //! [mesh Mesh Access]

  {
    //! [mesh Mesh Access 0]
    // Access number of cells
    unsigned num_entities = mesh_p->NumEntities(0);

    // Access number of edges
    unsigned num_nodes = mesh_p->NumEntities(1);

    //! [mesh Mesh Access 0]
  }

  {
    //! [mesh Entities]
    // Get entity pointer by index from a mesh
    const lf::mesh::Entity* entity = mesh_p->EntityByIndex(0, 0);

    // Iterate over all entities of co-dimension 0 (cells)
    for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
      // Do something with entity e.g.: print index of cell
      std::cout << mesh_p->Index(*cell) << std::endl;
    }
    //! [mesh Entities]

    //! [mesh Entities 0]
    // Get the co-dimension of the entity
    unsigned codim = entity->Codim();

    // Get the geometry of the entity, more details can be found
    // in the Geometry Quick reference
    const lf::geometry::Geometry* geometry = entity->Geometry();

    // Get the reference element of the entity (e.g. the unit triangle for a
    // general triangle)
    lf::base::RefEl ref_el = entity->RefEl();
    //! [mesh Entities 0]

    //! [mesh Entities 1]
    // Return all sub entities of this entity that have the given
    // co-dimension (w.r.t. this entity!)
    // For example, for a cell, the sub-entities of co-dimension 1
    std::span<const lf::mesh::Entity* const> sub_entities =
        entity->SubEntities(1);

    // If you require a vector of sub-entities, you can use the following code
    // snippet
    std::vector<const lf::mesh::Entity*> sub_entities_vec{
        entity->SubEntities(1).begin(), entity->SubEntities(1).end()};

    //! [mesh Entities 1]

    //! [mesh Entities 2]
    // return span of relative orientations of sub-entities of the next higher
    // co-dimension.
    std::span<const lf::mesh::Orientation> orientations =
        entity->RelativeOrientations();

    // If you require a vector of orientations, you can use the following code
    // snippet
    std::vector<lf::mesh::Orientation> orientations_vec{
        entity->RelativeOrientations().begin(),
        entity->RelativeOrientations().end()};
    //! [mesh Entities 2]
  }

  {
    //! [mesh Mesh Data Sets]
    // Create a mesh data set storing boolean values for each entity
    lf::mesh::utils::AllCodimMeshDataSet<bool> mesh_data_set =
        lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p);

    // Create a MDS storing boolean values for each entity of co-dimension 1
    // (edges)
    lf::mesh::utils::CodimMeshDataSet<bool> mesh_data_set_edges =
        lf::mesh::utils::CodimMeshDataSet<bool>(mesh_p, 1);

    // Mesh data sets can be initialized with a default value
    lf::mesh::utils::AllCodimMeshDataSet<bool> mesh_data_set_2 =
        lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p, false);

    // Access a (modifiable) value associated with an entity
    const lf::mesh::Entity* entity = mesh_p->EntityByIndex(0, 0);
    bool value = mesh_data_set(*entity);
    //! [mesh Mesh Data Sets]
  }

  {
    //! [mesh Mesh Data Sets 0]
    // Create MeshDataSet storing boolean values for each entity indicating if
    // the entity is on the boundary
    lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags{
        lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
    //! [mesh Mesh Data Sets 0]
  }

  {
    //! [mesh Mesh Functions]
    // Define a lambda function that takes a point in the reference element of
    // an entity and returns a value
    auto alpha = [](Eigen::Vector2d x) -> double { return 0.5 * x.norm(); };

    // Wrap the lambda function in a MeshFunctionGlobal object
    // Type of mesh function template parameter is automatically deduced
    lf::mesh::utils::MeshFunctionGlobal mesh_function =
        lf::mesh::utils::MeshFunctionGlobal(alpha);

    // Evaluate the mesh function on a cell
    const lf::mesh::Entity* cell = mesh_p->EntityByIndex(0, 0);
    std::vector<double> values =
        mesh_function(*cell, Eigen::Vector2d{0.5, 0.5});

    // We only asked for one point, so the vector has only one entry
    std::cout << "Value: " << values[0] << std::endl;  // Value: 0.790569
    //! [mesh Mesh Functions]
  }

  {
    //! [mesh Mesh Creation]
    // Generate a simple test mesh
    std::shared_ptr<lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
    //! [mesh Mesh Creation]
  }

  {
    //! [mesh Mesh Creation 0]
    // Read a mesh from a file
    std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader gmsh_reader =
        lf::io::GmshReader(std::move(mesh_factory), "triangle.msh");

    // Get pointer to the mesh
    std::shared_ptr<lf::mesh::Mesh> mesh = gmsh_reader.mesh();
    //! [mesh Mesh Creation 0]
  }

  {
    //! [mesh Mesh Refinement]
    std::shared_ptr<lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
    // Generate mesh hierarchy by uniform refinement with 6 levels

    std::shared_ptr<lf::refinement::MeshHierarchy> mesh_seq_p{
        lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 3)};

    // Access the mesh at level 3
    std::shared_ptr<lf::mesh::Mesh> mesh_level_3 = mesh_seq_p->getMesh(3);

    // We can refine a mesh further by calling
    mesh_seq_p->RefineRegular();

    // Access the mesh at level 4
    std::shared_ptr<lf::mesh::Mesh> mesh_level_4 = mesh_seq_p->getMesh(4);
    //! [mesh Mesh Refinement]
  }

  {
    //! [mesh Mesh Builder]
    // builder for a hybrid mesh in a world of dimension 2
    std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory =
        std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

    // Add points (needs to be done before adding entities)
    mesh_factory->AddPoint(Eigen::Vector2d{0, 0});    // (0)
    mesh_factory->AddPoint(Eigen::Vector2d{1, 0});    // (1)
    mesh_factory->AddPoint(Eigen::Vector2d{1, 1});    // (2)
    mesh_factory->AddPoint(Eigen::Vector2d{0, 1});    // (3)
    mesh_factory->AddPoint(Eigen::Vector2d{0.5, 1});  // (4)

    // Add a triangle
    // First set the coordinates of its nodes:
    Eigen::MatrixXd nodesOfTria(2, 3);
    nodesOfTria << 1, 1, 0.5, 0, 1, 1;
    mesh_factory->AddEntity(
        lf::base::RefEl::kTria(),  // we want a triangle
        std::array<lf::mesh::Mesh::size_type, 3>{
            {1, 2, 4}},  // indices of the nodes
        std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

    // Add a quadrilateral
    Eigen::MatrixXd nodesOfQuad(2, 4);
    nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
    mesh_factory->AddEntity(
        lf::base::RefEl::kQuad(),  // we want a quadrilateral
        std::array<lf::mesh::Mesh::size_type, 4>{
            {0, 1, 4, 3}},  // indices of the nodes
        std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));  // node coords

    // Build the mesh
    std::shared_ptr<lf::mesh::Mesh> mesh = mesh_factory->Build();
    //! [mesh Mesh Builder]
  }

  {
    //! [mesh Mesh Builder 0]
    // Create a LehrFEM++ square tensor product mesh
    lf::mesh::utils::TPTriagMeshBuilder builder(
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));

    // Set mesh parameters following the builder pattern
    // Domain is the unit square
    builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
        .setTopRightCorner(Eigen::Vector2d{99, 99})
        .setNumXCells(100)
        .setNumYCells(100);

    std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();
    //! [mesh Mesh Builder 0]
  }
}
