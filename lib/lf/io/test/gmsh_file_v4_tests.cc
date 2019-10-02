/**
 * @file
 * @brief Make sure that we can read a GmshFile Version 4
 * @author Raffael Casagrande
 * @date   2019-09-08 09:36:08
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>

namespace lf::io::test {

void checkPieceOfCake(const GMshFileV4& file) {
  EXPECT_EQ(file.version_number, "4.1");
  EXPECT_EQ(sizeof(std::size_t), file.size_t_size);
  EXPECT_EQ(file.physical_names.size(), 2);
  EXPECT_EQ(file.physical_names[0].dimension, 0);
  EXPECT_EQ(file.physical_names[0].physical_tag, 1);
  EXPECT_EQ(file.physical_names[0].name, "origin");
  EXPECT_EQ(file.physical_names[1].dimension, 1);
  EXPECT_EQ(file.physical_names[1].physical_tag, 2);
  EXPECT_EQ(file.physical_names[1].name, "arc");

  // entities:
  auto point_entities = std::get<0>(file.entities);
  auto curve_entities = std::get<1>(file.entities);
  auto surface_entities = std::get<2>(file.entities);
  auto volume_entities = std::get<3>(file.entities);

  EXPECT_EQ(point_entities.size(), 3);
  EXPECT_EQ(point_entities[0].tag, 1);
  EXPECT_EQ(point_entities[0].coord, Eigen::Vector3d::Zero());
  EXPECT_EQ(point_entities[0].physical_tags, std::vector<int>{1});
  EXPECT_EQ(point_entities[1].tag, 2);
  EXPECT_TRUE(point_entities[1].coord.isApprox(
      Eigen::Vector3d(1. / std::sqrt(2), 1. / std::sqrt(2), 0)));
  EXPECT_EQ(point_entities[1].physical_tags, std::vector<int>{});
  EXPECT_EQ(point_entities[2].tag, 3);
  EXPECT_TRUE(point_entities[2].coord.isApprox(Eigen::Vector3d(0, 1, 0)));
  EXPECT_EQ(point_entities[2].physical_tags, std::vector<int>{});

  EXPECT_EQ(curve_entities.size(), 3);
  EXPECT_EQ(curve_entities[0].tag, 1);
  EXPECT_TRUE(curve_entities[0].min_coord.isApprox(Eigen::Vector3d(
      -1.000000000028756e-007, -1.000000000028756e-007, -1e-007)));
  EXPECT_TRUE(curve_entities[0].max_coord.isApprox(
      Eigen::Vector3d(0.7071068811865474, 0.7071068811865474, 1e-007)));
  EXPECT_TRUE(curve_entities[0].physical_tags.empty());
  EXPECT_EQ(curve_entities[0].bounding_entities, (std::vector<int>{1, -2}));
  EXPECT_EQ(curve_entities[1].tag, 2);
  EXPECT_TRUE(curve_entities[1].min_coord.isApprox(
      Eigen::Vector3d(-9.999999983634211e-008, 0.7071066811865474, -1e-007)));
  EXPECT_TRUE(curve_entities[1].max_coord.isApprox(
      Eigen::Vector3d(0.7071068811865475, 1.0000001, 1e-007)));
  EXPECT_EQ(curve_entities[1].physical_tags, std::vector<int>{2});
  EXPECT_EQ(curve_entities[1].bounding_entities, (std::vector<int>{2, -3}));
  EXPECT_EQ(curve_entities[2].tag, 3);
  EXPECT_TRUE(curve_entities[2].min_coord.isApprox(
      Eigen::Vector3d(-1e-007, -9.999999994736442e-008, -1e-007)));
  EXPECT_TRUE(curve_entities[2].max_coord.isApprox(
      Eigen::Vector3d(1e-007, 1.0000001, 1e-007)));
  EXPECT_EQ(curve_entities[2].physical_tags, std::vector<int>{});
  EXPECT_EQ(curve_entities[2].bounding_entities, (std::vector<int>{3, -1}));

  EXPECT_EQ(surface_entities.size(), 1);
  EXPECT_EQ(surface_entities[0].tag, 1);
  EXPECT_TRUE(surface_entities[0].min_coord.isApprox(Eigen::Vector3d(
      -1.000000000028756e-007, -9.999999994736442e-008, -1e-007)));
  EXPECT_TRUE(surface_entities[0].max_coord.isApprox(
      Eigen::Vector3d(0.7071068811865475, 1.0000001, 1e-007)));
  EXPECT_EQ(surface_entities[0].physical_tags, std::vector<int>{3});
  EXPECT_EQ(surface_entities[0].bounding_entities, (std::vector<int>{1, 2, 3}));
  EXPECT_EQ(volume_entities.size(), 0);

  // partitioned entities:
  auto pe = file.partitioned_entities;
  EXPECT_EQ(pe.num_partitions, 2);
  EXPECT_EQ(pe.ghost_entities.size(), 2);
  EXPECT_EQ(pe.ghost_entities[0].tag, 4);
  EXPECT_EQ(pe.ghost_entities[0].partition, 1);
  EXPECT_EQ(pe.ghost_entities[1].tag, 5);
  EXPECT_EQ(pe.ghost_entities[1].partition, 2);

  auto ppe = std::get<0>(pe.partitioned_entities);
  auto pce = std::get<1>(pe.partitioned_entities);
  auto pse = std::get<2>(pe.partitioned_entities);
  auto pve = std::get<3>(pe.partitioned_entities);
  EXPECT_EQ(ppe.size(), 5);
  EXPECT_EQ(ppe[0].tag, 4);
  EXPECT_EQ(ppe[0].parent_dim, 0);
  EXPECT_EQ(ppe[0].parent_tag, 1);
  EXPECT_EQ(ppe[0].partitions, std::vector<int>{2});
  EXPECT_TRUE(ppe[0].coord.isZero());
  EXPECT_EQ(ppe[0].physical_tags, std::vector<int>{1});
  EXPECT_EQ(ppe[1].tag, 5);
  EXPECT_EQ(ppe[1].parent_dim, 0);
  EXPECT_EQ(ppe[1].parent_tag, 2);
  EXPECT_EQ(ppe[1].partitions, std::vector<int>{2});
  EXPECT_TRUE(ppe[1].coord.isApprox(
      Eigen::Vector3d(0.7071067811865475, 0.7071067811865475, 0)));
  EXPECT_EQ(ppe[1].physical_tags, std::vector<int>{});

  EXPECT_EQ(pce.size(), 5);
  EXPECT_EQ(pce[0].tag, 4);
  EXPECT_EQ(pce[0].parent_dim, 1);
  EXPECT_EQ(pce[0].parent_tag, 1);
  EXPECT_EQ(pce[0].partitions, std::vector<int>{2});
  EXPECT_TRUE(pce[0].min_coord.isZero());
  EXPECT_TRUE(pce[0].max_coord.isApprox(
      Eigen::Vector3d(0.7071067811865475, 0.7071067811865475, 0)));
  EXPECT_TRUE(pce[0].physical_tags.empty());
  EXPECT_EQ(pce[0].bounding_entities, (std::vector<int>{8, -5}));
  EXPECT_EQ(pce[1].tag, 5);
  EXPECT_EQ(pce[1].parent_dim, 1);
  EXPECT_EQ(pce[1].parent_tag, 2);
  EXPECT_EQ(pce[1].partitions, std::vector<int>{2});
  EXPECT_TRUE(pce[1].min_coord.isApprox(
      Eigen::Vector3d(0.3826834323650894, 0.7071067811865475, 0)));
  EXPECT_TRUE(pce[1].max_coord.isApprox(
      Eigen::Vector3d(0.7071067811865475, 0.9238795325112867, 0)));
  EXPECT_EQ(pce[1].physical_tags, std::vector<int>{2});
  EXPECT_EQ(pce[1].bounding_entities, (std::vector<int>{5, -7}));

  EXPECT_EQ(pse.size(), 2);
  EXPECT_EQ(pse[0].tag, 2);
  EXPECT_EQ(pse[0].parent_dim, 2);
  EXPECT_EQ(pse[0].parent_tag, 1);
  EXPECT_EQ(pse[0].partitions, std::vector<int>{2});
  EXPECT_TRUE(pse[0].min_coord.isZero());
  EXPECT_TRUE(pse[0].max_coord.isApprox(
      Eigen::Vector3d(0.7071067811865475, 0.9238795325112867, 0)));
  EXPECT_EQ(pse[0].physical_tags, std::vector<int>{3});
  // we have to use is_permutation because binary format permutes them
  // differently than text format
  EXPECT_TRUE(std::is_permutation(pse[0].bounding_entities.begin(),
                                  pse[0].bounding_entities.end(),
                                  (std::vector<int>{5, 4, -8}).begin()));
  EXPECT_EQ(pse[1].tag, 3);
  EXPECT_EQ(pse[1].parent_dim, 2);
  EXPECT_EQ(pse[1].parent_tag, 1);
  EXPECT_EQ(pse[1].partitions, std::vector<int>{1});
  EXPECT_TRUE(pse[1].min_coord.isZero());
  EXPECT_TRUE(
      pse[1].max_coord.isApprox(Eigen::Vector3d(0.3826834323650894, 1, 0)));
  EXPECT_EQ(pse[1].physical_tags, std::vector<int>{3});
  // we have to use is_permutation because binary format permutes them
  // differently than text format
  EXPECT_TRUE(std::is_permutation(pse[1].bounding_entities.begin(),
                                  pse[1].bounding_entities.end(),
                                  (std::vector<int>{6, 7, 8}).begin()));

  // nodes
  EXPECT_EQ(file.nodes.num_nodes, 4);
  EXPECT_EQ(file.nodes.min_node_tag, 1);
  EXPECT_EQ(file.nodes.max_node_tag, 4);
  auto nb = file.nodes.node_blocks;
  EXPECT_EQ(nb.size(), 12);
  EXPECT_EQ(nb[0].entity_dim, 0);
  EXPECT_EQ(nb[0].entity_tag, 4);
  EXPECT_FALSE(nb[0].parametric);
  EXPECT_EQ(nb[0].nodes.size(), 1);
  EXPECT_EQ(nb[0].nodes[0].first, 1);
  EXPECT_TRUE(nb[0].nodes[0].second.isZero());
  EXPECT_EQ(nb[1].entity_dim, 0);
  EXPECT_EQ(nb[1].entity_tag, 5);
  EXPECT_FALSE(nb[1].parametric);
  EXPECT_EQ(nb[1].nodes.size(), 1);
  EXPECT_EQ(nb[1].nodes[0].first, 2);
  EXPECT_TRUE(nb[1].nodes[0].second.isApprox(
      Eigen::Vector3d(0.7071067811865475, 0.7071067811865475, 0)));
  EXPECT_EQ(nb[2].entity_dim, 0);
  EXPECT_EQ(nb[2].entity_tag, 6);
  EXPECT_FALSE(nb[2].parametric);
  EXPECT_EQ(nb[2].nodes.size(), 1);
  EXPECT_EQ(nb[2].nodes[0].first, 3);
  EXPECT_TRUE(nb[2].nodes[0].second.isApprox(Eigen::Vector3d(0, 1, 0)));
  EXPECT_EQ(nb[3].entity_dim, 0);
  EXPECT_EQ(nb[3].entity_tag, 7);
  EXPECT_FALSE(nb[3].parametric);
  EXPECT_EQ(nb[3].nodes.size(), 1);
  EXPECT_EQ(nb[3].nodes[0].first, 4);
  EXPECT_TRUE(nb[3].nodes[0].second.isApprox(
      Eigen::Vector3d(0.3826834323650894, 0.9238795325112867, 0)));
  EXPECT_EQ(nb[4].entity_dim, 0);
  EXPECT_EQ(nb[4].entity_tag, 8);
  EXPECT_FALSE(nb[4].parametric);
  EXPECT_EQ(nb[4].nodes.size(), 0);
  EXPECT_EQ(nb[5].entity_dim, 1);
  EXPECT_EQ(nb[5].entity_tag, 4);
  EXPECT_FALSE(nb[5].parametric);
  EXPECT_EQ(nb[5].nodes.size(), 0);
  EXPECT_EQ(nb[6].entity_dim, 1);
  EXPECT_EQ(nb[6].entity_tag, 5);
  EXPECT_FALSE(nb[6].parametric);
  EXPECT_EQ(nb[6].nodes.size(), 0);
  EXPECT_EQ(nb[7].entity_dim, 1);
  EXPECT_EQ(nb[7].entity_tag, 6);
  EXPECT_FALSE(nb[7].parametric);
  EXPECT_EQ(nb[7].nodes.size(), 0);
  EXPECT_EQ(nb[8].entity_dim, 1);
  EXPECT_EQ(nb[8].entity_tag, 7);
  EXPECT_FALSE(nb[8].parametric);
  EXPECT_EQ(nb[8].nodes.size(), 0);
  EXPECT_EQ(nb[9].entity_dim, 1);
  EXPECT_EQ(nb[9].entity_tag, 8);
  EXPECT_FALSE(nb[9].parametric);
  EXPECT_EQ(nb[9].nodes.size(), 0);
  EXPECT_EQ(nb[10].entity_dim, 2);
  EXPECT_EQ(nb[10].entity_tag, 2);
  EXPECT_FALSE(nb[10].parametric);
  EXPECT_EQ(nb[10].nodes.size(), 0);
  EXPECT_EQ(nb[11].entity_dim, 2);
  EXPECT_EQ(nb[11].entity_tag, 3);
  EXPECT_FALSE(nb[11].parametric);
  EXPECT_EQ(nb[11].nodes.size(), 0);

  // elements:
  EXPECT_EQ(file.elements.num_elements, 7);
  EXPECT_EQ(file.elements.min_element_tag, 1);
  EXPECT_EQ(file.elements.max_element_tag, 11);
  auto& eb = file.elements.element_blocks;
  EXPECT_EQ(eb.size(), 7);
  EXPECT_EQ(eb[0].dimension, 0);
  EXPECT_EQ(eb[0].entity_tag, 4);
  EXPECT_EQ(eb[0].element_type, GMshFileV4::ElementType::POINT);
  EXPECT_EQ(eb[0].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[0].elements[0]), 1);
  EXPECT_EQ(std::get<1>(eb[0].elements[0]), std::vector<std::size_t>{1});
  EXPECT_EQ(eb[1].dimension, 0);
  EXPECT_EQ(eb[1].entity_tag, 7);
  EXPECT_EQ(eb[1].element_type, GMshFileV4::ElementType::POINT);
  EXPECT_EQ(eb[1].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[1].elements[0]), 11);
  EXPECT_EQ(std::get<1>(eb[1].elements[0]), std::vector<std::size_t>{4});
  EXPECT_EQ(eb[2].dimension, 1);
  EXPECT_EQ(eb[2].entity_tag, 5);
  EXPECT_EQ(eb[2].element_type, GMshFileV4::ElementType::EDGE2);
  EXPECT_EQ(eb[2].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[2].elements[0]), 2);
  EXPECT_EQ(std::get<1>(eb[2].elements[0]), (std::vector<std::size_t>{2, 4}));
  EXPECT_EQ(eb[3].dimension, 1);
  EXPECT_EQ(eb[3].entity_tag, 6);
  EXPECT_EQ(eb[3].element_type, GMshFileV4::ElementType::EDGE2);
  EXPECT_EQ(eb[3].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[3].elements[0]), 3);
  EXPECT_EQ(std::get<1>(eb[3].elements[0]), (std::vector<std::size_t>{4, 3}));
  EXPECT_EQ(eb[4].dimension, 1);
  EXPECT_EQ(eb[4].entity_tag, 8);
  EXPECT_EQ(eb[4].element_type, GMshFileV4::ElementType::EDGE2);
  EXPECT_EQ(eb[4].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[4].elements[0]), 10);
  EXPECT_EQ(std::get<1>(eb[4].elements[0]), (std::vector<std::size_t>{1, 4}));
  EXPECT_EQ(eb[5].dimension, 2);
  EXPECT_EQ(eb[5].entity_tag, 2);
  EXPECT_EQ(eb[5].element_type, GMshFileV4::ElementType::TRIA3);
  EXPECT_EQ(eb[5].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[5].elements[0]), 4);
  EXPECT_EQ(std::get<1>(eb[5].elements[0]),
            (std::vector<std::size_t>{1, 2, 4}));
  EXPECT_EQ(eb[6].dimension, 2);
  EXPECT_EQ(eb[6].entity_tag, 3);
  EXPECT_EQ(eb[6].element_type, GMshFileV4::ElementType::TRIA3);
  EXPECT_EQ(eb[6].elements.size(), 1);
  EXPECT_EQ(std::get<0>(eb[6].elements[0]), 5);
  EXPECT_EQ(std::get<1>(eb[6].elements[0]),
            (std::vector<std::size_t>{1, 4, 3}));

  // periodic
  Eigen::Matrix4d affine_transform;
  affine_transform << 0.7071067811865476, 0.7071067811865475, 0, 0,
      -0.7071067811865475, 0.7071067811865476, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
  EXPECT_EQ(file.periodic_links.size(), 2);
  EXPECT_EQ(file.periodic_links[0].dimension, 0);
  EXPECT_EQ(file.periodic_links[0].entity_tag_slave, 2);
  EXPECT_EQ(file.periodic_links[0].entity_tag_master, 3);
  EXPECT_TRUE(file.periodic_links[0].affine_transform);
  EXPECT_TRUE(
      file.periodic_links[0].affine_transform->isApprox(affine_transform));
  EXPECT_EQ(file.periodic_links[0].node_mapping.size(), 1);
  EXPECT_EQ(file.periodic_links[0].node_mapping[0].first, 2);
  EXPECT_EQ(file.periodic_links[0].node_mapping[0].second, 3);
  EXPECT_EQ(file.periodic_links[1].dimension, 1);
  EXPECT_EQ(file.periodic_links[1].entity_tag_slave, 1);
  EXPECT_EQ(file.periodic_links[1].entity_tag_master, 3);
  EXPECT_TRUE(file.periodic_links[1].affine_transform);
  EXPECT_TRUE(
      file.periodic_links[1].affine_transform->isApprox(affine_transform));
  EXPECT_EQ(file.periodic_links[1].node_mapping.size(), 2);
  EXPECT_EQ(file.periodic_links[1].node_mapping[0].first, 1);
  EXPECT_EQ(file.periodic_links[1].node_mapping[0].second, 1);
  EXPECT_EQ(file.periodic_links[1].node_mapping[1].first, 2);
  EXPECT_EQ(file.periodic_links[1].node_mapping[1].second, 3);

  // ghost elements:
  EXPECT_EQ(file.ghost_elements.size(), 2);
  if (file.ghost_elements[0].element_tag == 5) {
    EXPECT_EQ(file.ghost_elements[0].element_tag, 5);
    EXPECT_EQ(file.ghost_elements[0].partition_tag, 1);
    EXPECT_EQ(file.ghost_elements[0].ghost_partition_tags, std::vector<int>{2});
    EXPECT_EQ(file.ghost_elements[1].element_tag, 4);
    EXPECT_EQ(file.ghost_elements[1].partition_tag, 2);
    EXPECT_EQ(file.ghost_elements[1].ghost_partition_tags, std::vector<int>{1});
  } else {
    EXPECT_EQ(file.ghost_elements[1].element_tag, 5);
    EXPECT_EQ(file.ghost_elements[1].partition_tag, 1);
    EXPECT_EQ(file.ghost_elements[1].ghost_partition_tags, std::vector<int>{2});
    EXPECT_EQ(file.ghost_elements[0].element_tag, 4);
    EXPECT_EQ(file.ghost_elements[0].partition_tag, 2);
    EXPECT_EQ(file.ghost_elements[0].ghost_partition_tags, std::vector<int>{1});
  }
}

void checkTwoElementHybrid(GMshFileV4 file) {
  // Header
  EXPECT_EQ(file.version_number, "4.1");
  EXPECT_EQ(file.size_t_size, sizeof(std::size_t));

  // Physical names:
  EXPECT_EQ(file.physical_names.size(), 6);
  EXPECT_EQ(file.physical_names[0].dimension, 0);
  EXPECT_EQ(file.physical_names[0].physical_tag, 1);
  EXPECT_EQ(file.physical_names[0].name, "physicalEntity1");
  EXPECT_EQ(file.physical_names[1].dimension, 0);
  EXPECT_EQ(file.physical_names[1].physical_tag, 2);
  EXPECT_EQ(file.physical_names[1].name, "physicalEntity2");
  EXPECT_EQ(file.physical_names[2].dimension, 1);
  EXPECT_EQ(file.physical_names[2].physical_tag, 4);
  EXPECT_EQ(file.physical_names[2].name, "diagonal");
  EXPECT_EQ(file.physical_names[3].dimension, 2);
  EXPECT_EQ(file.physical_names[3].physical_tag, 1);
  EXPECT_EQ(file.physical_names[3].name, "physicalEntity1");
  EXPECT_EQ(file.physical_names[4].dimension, 2);
  EXPECT_EQ(file.physical_names[4].physical_tag, 3);
  EXPECT_EQ(file.physical_names[4].name, "physicalEntity3");
  EXPECT_EQ(file.physical_names[5].dimension, 2);
  EXPECT_EQ(file.physical_names[5].physical_tag, 5);
  EXPECT_EQ(file.physical_names[5].name, "square");

  // Entities:
  auto& points = std::get<0>(file.entities);
  EXPECT_EQ(points.size(), 5);
  EXPECT_EQ(points[0].tag, 1);
  EXPECT_EQ(points[0].coord, Eigen::Vector3d(0, 0, 0));
  EXPECT_EQ(points[0].physical_tags, (std::vector<int>{1, 2}));
  EXPECT_EQ(points[1].tag, 2);
  EXPECT_EQ(points[1].coord, Eigen::Vector3d(1, 0, 0));
  EXPECT_TRUE(points[1].physical_tags.empty());
  EXPECT_EQ(points[2].tag, 3);
  EXPECT_EQ(points[2].coord, Eigen::Vector3d(2, 0, 0));
  EXPECT_TRUE(points[2].physical_tags.empty());
  EXPECT_EQ(points[3].tag, 4);
  EXPECT_EQ(points[3].coord, Eigen::Vector3d(0, 1, 0));
  EXPECT_TRUE(points[3].physical_tags.empty());
  EXPECT_EQ(points[4].tag, 5);
  EXPECT_EQ(points[4].coord, Eigen::Vector3d(1, 1, 0));
  EXPECT_TRUE(points[4].physical_tags.empty());

  auto& curves = std::get<1>(file.entities);
  EXPECT_EQ(curves.size(), 6);
  EXPECT_EQ(curves[0].tag, 1);
  EXPECT_TRUE(curves[0].min_coord.isApprox(
      Eigen::Vector3d{-9.999999994736442e-008, -1e-007, -1e-007}));
  EXPECT_TRUE(
      curves[0].max_coord.isApprox(Eigen::Vector3d(1.0000001, 1e-007, 1e-007)));
  EXPECT_TRUE(curves[0].physical_tags.empty());
  EXPECT_EQ(curves[0].bounding_entities, (std::vector<int>{1, -2}));
  EXPECT_EQ(curves[1].tag, 2);
  EXPECT_TRUE(curves[1].min_coord.isApprox(
      Eigen::Vector3d{0.9999999000000001, -1e-007, -1e-007}));
  EXPECT_TRUE(
      curves[1].max_coord.isApprox(Eigen::Vector3d(2.0000001, 1e-007, 1e-007)));
  EXPECT_TRUE(curves[1].physical_tags.empty());
  EXPECT_EQ(curves[1].bounding_entities, (std::vector<int>{2, -3}));
  EXPECT_EQ(curves[2].tag, 3);
  EXPECT_TRUE(curves[2].min_coord.isApprox(
      Eigen::Vector3d{-1e-007, -9.999999994736442e-008, -1e-007}));
  EXPECT_TRUE(
      curves[2].max_coord.isApprox(Eigen::Vector3d(1e-007, 1.0000001, 1e-007)));
  EXPECT_TRUE(curves[2].physical_tags.empty());
  EXPECT_EQ(curves[2].bounding_entities, (std::vector<int>{1, -4}));
  EXPECT_EQ(curves[3].tag, 4);
  EXPECT_TRUE(curves[3].min_coord.isApprox(
      Eigen::Vector3d{0.9999999000000001, -9.999999994736442e-008, -1e-007}));
  EXPECT_TRUE(curves[3].max_coord.isApprox(
      Eigen::Vector3d(1.0000001, 1.0000001, 1e-007)));
  EXPECT_TRUE(curves[3].physical_tags.empty());
  EXPECT_EQ(curves[3].bounding_entities, (std::vector<int>{2, -5}));
  EXPECT_EQ(curves[4].tag, 5);
  EXPECT_TRUE(curves[4].min_coord.isApprox(
      Eigen::Vector3d{-9.999999994736442e-008, 0.9999999000000001, -1e-007}));
  EXPECT_TRUE(curves[4].max_coord.isApprox(
      Eigen::Vector3d(1.0000001, 1.0000001, 1e-007)));
  EXPECT_TRUE(curves[4].physical_tags.empty());
  EXPECT_EQ(curves[4].bounding_entities, (std::vector<int>{4, -5}));
  EXPECT_EQ(curves[5].tag, 6);
  EXPECT_TRUE(curves[5].min_coord.isApprox(
      Eigen::Vector3d{0.9999999000000001, -9.999999994736442e-008, -1e-007}));
  EXPECT_TRUE(curves[5].max_coord.isApprox(
      Eigen::Vector3d(2.0000001, 1.0000001, 1e-007)));
  EXPECT_EQ(curves[5].physical_tags, std::vector<int>{4});
  EXPECT_EQ(curves[5].bounding_entities, (std::vector<int>{3, -5}));

  auto& surfaces = std::get<2>(file.entities);
  EXPECT_EQ(surfaces.size(), 2);
  EXPECT_EQ(surfaces[0].tag, 1);
  EXPECT_TRUE(surfaces[0].min_coord.isApprox(Eigen::Vector3d(
      -9.999999994736442e-008, -9.999999994736442e-008, -1e-007)));
  EXPECT_TRUE(surfaces[0].max_coord.isApprox(
      Eigen::Vector3d(1.0000001, 1.0000001, 1e-007)));
  EXPECT_EQ(surfaces[0].physical_tags, std::vector<int>{5});
  EXPECT_EQ(surfaces[0].bounding_entities, (std::vector<int>{3, 5, -4, -1}));
  EXPECT_EQ(surfaces[1].tag, 2);
  EXPECT_TRUE(surfaces[1].min_coord.isApprox(
      Eigen::Vector3d(0.9999999000000001, -9.999999994736442e-008, -1e-007)));
  EXPECT_TRUE(surfaces[1].max_coord.isApprox(
      Eigen::Vector3d(2.0000001, 1.0000001, 1e-007)));
  EXPECT_EQ(surfaces[1].physical_tags, (std::vector<int>{1, 3}));
  EXPECT_EQ(surfaces[1].bounding_entities, (std::vector<int>{4, -6, -2}));
  EXPECT_TRUE(std::get<3>(file.entities).empty());

  // partitioned entities
  EXPECT_EQ(file.partitioned_entities.num_partitions, 0);
  EXPECT_TRUE(file.partitioned_entities.ghost_entities.empty());
  EXPECT_TRUE(
      std::get<0>(file.partitioned_entities.partitioned_entities).empty());
  EXPECT_TRUE(
      std::get<1>(file.partitioned_entities.partitioned_entities).empty());
  EXPECT_TRUE(
      std::get<2>(file.partitioned_entities.partitioned_entities).empty());
  EXPECT_TRUE(
      std::get<3>(file.partitioned_entities.partitioned_entities).empty());

  // nodes
  using vnode_t = std::vector<std::pair<std::size_t, Eigen::Vector3d>>;
  EXPECT_EQ(file.nodes.num_nodes, 5);
  EXPECT_EQ(file.nodes.min_node_tag, 1);
  EXPECT_EQ(file.nodes.max_node_tag, 5);
  EXPECT_EQ(file.nodes.node_blocks[0].entity_dim, 0);
  EXPECT_EQ(file.nodes.node_blocks[0].entity_tag, 1);
  EXPECT_FALSE(file.nodes.node_blocks[0].parametric);
  EXPECT_EQ(file.nodes.node_blocks[0].nodes,
            (vnode_t{{1, Eigen::Vector3d(0, 0, 0)}}));
  EXPECT_EQ(file.nodes.node_blocks[1].entity_dim, 0);
  EXPECT_EQ(file.nodes.node_blocks[1].entity_tag, 2);
  EXPECT_FALSE(file.nodes.node_blocks[1].parametric);
  EXPECT_EQ(file.nodes.node_blocks[1].nodes,
            (vnode_t{{2, Eigen::Vector3d(1, 0, 0)}}));
  EXPECT_EQ(file.nodes.node_blocks[2].entity_dim, 0);
  EXPECT_EQ(file.nodes.node_blocks[2].entity_tag, 3);
  EXPECT_FALSE(file.nodes.node_blocks[2].parametric);
  EXPECT_EQ(file.nodes.node_blocks[2].nodes,
            (vnode_t{{3, Eigen::Vector3d(2, 0, 0)}}));
  EXPECT_EQ(file.nodes.node_blocks[3].entity_dim, 0);
  EXPECT_EQ(file.nodes.node_blocks[3].entity_tag, 4);
  EXPECT_FALSE(file.nodes.node_blocks[3].parametric);
  EXPECT_EQ(file.nodes.node_blocks[3].nodes,
            (vnode_t{{4, Eigen::Vector3d(0, 1, 0)}}));
  EXPECT_EQ(file.nodes.node_blocks[4].entity_dim, 0);
  EXPECT_EQ(file.nodes.node_blocks[4].entity_tag, 5);
  EXPECT_FALSE(file.nodes.node_blocks[4].parametric);
  EXPECT_EQ(file.nodes.node_blocks[4].nodes,
            (vnode_t{{5, Eigen::Vector3d(1, 1, 0)}}));
  EXPECT_EQ(file.nodes.node_blocks[5].entity_dim, 1);
  EXPECT_EQ(file.nodes.node_blocks[5].entity_tag, 6);
  EXPECT_FALSE(file.nodes.node_blocks[5].parametric);
  EXPECT_TRUE(file.nodes.node_blocks[5].nodes.empty());
  EXPECT_EQ(file.nodes.node_blocks[6].entity_dim, 2);
  EXPECT_EQ(file.nodes.node_blocks[6].entity_tag, 1);
  EXPECT_FALSE(file.nodes.node_blocks[6].parametric);
  EXPECT_TRUE(file.nodes.node_blocks[6].nodes.empty());
  EXPECT_EQ(file.nodes.node_blocks[7].entity_dim, 2);
  EXPECT_EQ(file.nodes.node_blocks[7].entity_tag, 2);
  EXPECT_FALSE(file.nodes.node_blocks[7].parametric);
  EXPECT_TRUE(file.nodes.node_blocks[7].nodes.empty());

  // elements
  using element_t =
      std::vector<std::tuple<std::size_t, std::vector<std::size_t>>>;
  EXPECT_EQ(file.elements.num_elements, 4);
  EXPECT_EQ(file.elements.min_element_tag, 1);
  EXPECT_EQ(file.elements.max_element_tag, 4);
  EXPECT_EQ(file.elements.element_blocks.size(), 4);
  EXPECT_EQ(file.elements.element_blocks[0].dimension, 0);
  EXPECT_EQ(file.elements.element_blocks[0].entity_tag, 1);
  EXPECT_EQ(file.elements.element_blocks[0].element_type,
            GMshFileV4::ElementType::POINT);
  EXPECT_EQ(file.elements.element_blocks[0].elements, (element_t{{1, {1}}}));
  EXPECT_EQ(file.elements.element_blocks[1].dimension, 1);
  EXPECT_EQ(file.elements.element_blocks[1].entity_tag, 6);
  EXPECT_EQ(file.elements.element_blocks[1].element_type,
            GMshFileV4::ElementType::EDGE2);
  EXPECT_EQ(file.elements.element_blocks[1].elements, (element_t{{2, {3, 5}}}));
  EXPECT_EQ(file.elements.element_blocks[2].dimension, 2);
  EXPECT_EQ(file.elements.element_blocks[2].entity_tag, 1);
  EXPECT_EQ(file.elements.element_blocks[2].element_type,
            GMshFileV4::ElementType::QUAD4);
  EXPECT_EQ(file.elements.element_blocks[2].elements,
            (element_t{{3, {1, 4, 5, 2}}}));
  EXPECT_EQ(file.elements.element_blocks[3].dimension, 2);
  EXPECT_EQ(file.elements.element_blocks[3].entity_tag, 2);
  EXPECT_EQ(file.elements.element_blocks[3].element_type,
            GMshFileV4::ElementType::TRIA3);
  EXPECT_EQ(file.elements.element_blocks[3].elements,
            (element_t{{4, {3, 5, 2}}}));
}

TEST(lf_io_gmsh_file_v4, readPieceOfCake) {
  auto filename = test_utils::getMeshPath("piece_of_cake.msh");
  auto mshfile = ReadGmshFile(filename);
  EXPECT_TRUE(std::holds_alternative<GMshFileV4>(mshfile));
  auto filev4 = std::get<GMshFileV4>(mshfile);
  EXPECT_FALSE(filev4.is_binary);
  checkPieceOfCake(filev4);
}

TEST(lf_io_gmsh_file_v4, readPieceOfCakeBinary) {
  auto filename = test_utils::getMeshPath("piece_of_cake_binary.msh");
  auto mshfile = ReadGmshFile(filename);
  EXPECT_TRUE(std::holds_alternative<GMshFileV4>(mshfile));
  auto filev4 = std::get<GMshFileV4>(mshfile);
  EXPECT_TRUE(filev4.is_binary);
  checkPieceOfCake(filev4);
}

TEST(lf_io_gmsh_file_v4, readTwoElementHybrid2d) {
  auto filename = test_utils::getMeshPath("two_element_hybrid_2d_v4.msh");
  auto mshfile = ReadGmshFile(filename);
  EXPECT_TRUE(std::holds_alternative<GMshFileV4>(mshfile));
  auto filev4 = std::get<GMshFileV4>(mshfile);
  EXPECT_FALSE(filev4.is_binary);
  checkTwoElementHybrid(filev4);
}

TEST(lf_io_gmsh_file_v4, readTwoElementHybrid2dBinary) {
  auto filename =
      test_utils::getMeshPath("two_element_hybrid_2d_v4_binary.msh");
  auto mshfile = ReadGmshFile(filename);
  EXPECT_TRUE(std::holds_alternative<GMshFileV4>(mshfile));
  auto filev4 = std::get<GMshFileV4>(mshfile);
  EXPECT_TRUE(filev4.is_binary);
  checkTwoElementHybrid(filev4);
}

}  // namespace lf::io::test
