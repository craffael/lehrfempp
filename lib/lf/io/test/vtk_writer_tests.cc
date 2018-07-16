/**
 * @file
 * @brief Test the vtk writer implementation.
 * @author Raffael Casagrande
 * @date   2018-07-14 07:52:25
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>

namespace lf::io::test {

template <class T>
using FieldDataArray = VtkFile::FieldDataArray<T>;

template <class T>
using ScalarData = VtkFile::ScalarData<T>;

template <class T>
using VectorData = VtkFile::VectorData<T>;

TEST(lf_io_VtkWriter, writeVtkFile) {
  VtkFile vtk_file;
  vtk_file.header = "this is my test header :)";
  vtk_file.unstructured_grid.points = {
      {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  vtk_file.unstructured_grid.cells = {{0, 1, 3, 4}, {1, 2, 3}};
  vtk_file.unstructured_grid.cell_types = {VtkFile::CellType::VTK_QUAD,
                                           VtkFile::CellType::VTK_TRIANGLE};

  vtk_file.field_data.push_back(FieldDataArray<int>("array0", {0, 1, 2}));
  vtk_file.field_data.push_back(FieldDataArray<float>("array1", {-1, 0, 1}));
  vtk_file.field_data.push_back(FieldDataArray<double>("array2", {0, 1, 2}));

  vtk_file.point_data.push_back(ScalarData<char>("char", {-2, -1, 0, 1, 2}));
  vtk_file.point_data.push_back(
      ScalarData<unsigned char>("uchar", {1, 2, 3, 4, 5}));
  vtk_file.point_data.push_back(ScalarData<short>("short", {-2, -1, 0, 1, 2}));
  vtk_file.point_data.push_back(
      ScalarData<unsigned short>("ushort", {1, 2, 3, 4, 5}));
  vtk_file.point_data.push_back(ScalarData<int>("int", {-2, -1, 0, 1, 2}));
  vtk_file.point_data.push_back(
      ScalarData<unsigned int>("uint", {1, 2, 3, 4, 5}));
  // vtk_file.point_data.push_back(ScalarData<long>("long", {-2, -1, 0, 1, 2}));
  // vtk_file.point_data.push_back(
  //     ScalarData<unsigned long>("ulong", {1, 2, 3, 4, 5}));
  vtk_file.point_data.push_back(ScalarData<float>("float", {-2, -1, 0, 1, 2}));
  vtk_file.point_data.push_back(
      ScalarData<double>("double", {-2, -1, 0, 1, 2}));

  vtk_file.point_data.push_back(VectorData<float>(
      "vfloat", {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {4, 0, 0}}));
  vtk_file.point_data.push_back(VectorData<double>(
      "vdouble", {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}, {4, 0, 0}}));

  vtk_file.cell_data.push_back(ScalarData<char>("char", {0, 1}));
  vtk_file.cell_data.push_back(ScalarData<unsigned char>("uchar", {0, 1}));
  vtk_file.cell_data.push_back(ScalarData<short>("short", {0, 1}));
  vtk_file.cell_data.push_back(ScalarData<unsigned short>("ushort", {0, 1}));
  vtk_file.cell_data.push_back(ScalarData<int>("int", {0, 1}));
  vtk_file.cell_data.push_back(ScalarData<unsigned int>("uint", {0, 1}));
  // vtk_file.cell_data.push_back(ScalarData<long>("long", {0, 1}));
  // vtk_file.cell_data.push_back(ScalarData<unsigned long>("ulong", {0, 1}));
  vtk_file.cell_data.push_back(
      VectorData<float>("vfloat", {{0, 0, 0}, {-1, 0, 0}}));
  vtk_file.cell_data.push_back(
      VectorData<double>("vdouble", {{0, 0, 0}, {-1, 0, 0}}));

  WriteToFile(vtk_file, "all_features.vtk");

  vtk_file.format = VtkFile::Format::BINARY;
  WriteToFile(vtk_file, "all_features_binary.vtk");
}

TEST(lf_io_VtkWriter, onlyMesh) {
  VtkFile vtk_file;
  vtk_file.header = "this is my test header :)";
  vtk_file.format = VtkFile::Format::ASCII;
  vtk_file.unstructured_grid.points = {
      {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {0, 1, 0}};
  vtk_file.unstructured_grid.cells = {{0, 1, 3, 4}, {1, 2, 3}};
  vtk_file.unstructured_grid.cell_types = {VtkFile::CellType::VTK_QUAD,
                                           VtkFile::CellType::VTK_TRIANGLE};

  WriteToFile(vtk_file, "only_mesh.vtk");

  vtk_file.format = VtkFile::Format::BINARY;
  WriteToFile(vtk_file, "only_mesh_binary.vtk");
}

}  // namespace lf::io::test
