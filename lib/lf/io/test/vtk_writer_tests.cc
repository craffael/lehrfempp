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

 TEST(lf_io_VtkWriter, writeVtkFile) {
   VtkFile vtk_file;
   vtk_file.header = "this is my test header :)";
   vtk_file.unstructured_grid.points = {
       {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {0, 1, 0}};
   vtk_file.unstructured_grid.cells = {{0, 1, 3, 4}, {1, 2, 3}};
   vtk_file.unstructured_grid.cell_types = {VtkFile::CellType::VTK_QUAD,
                                            VtkFile::CellType::VTK_TRIANGLE};

   vtk_file.point_data.scalar_int_data.resize(1);
   vtk_file.point_data.scalar_int_data[0].data_name = "continuous_sequence";
   vtk_file.point_data.scalar_int_data[0].data = {1, 2, 3, 4, 5};

   vtk_file.point_data.scalar_unsigned_int_data.resize(2);
   vtk_file.point_data.scalar_unsigned_int_data[0].data_name = "all1";
   vtk_file.point_data.scalar_unsigned_int_data[0].data = {1, 1, 1, 1, 1};
   vtk_file.point_data.scalar_unsigned_int_data[1].data_name = "all2";
   vtk_file.point_data.scalar_unsigned_int_data[1].data = {2, 2, 2, 2, 2};

   WriteToFile(vtk_file, "test.vtk");
 }

 // TEST(lf_io_VtkWriter, onlyMesh) {
 //   VtkFile vtk_file;
 //   vtk_file.header = "this is my test header :)";
 //   vtk_file.format = VtkFile::Format::BINARY;
 //   vtk_file.unstructured_grid.points = {
 //       {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {0, 1, 0}};
 //   vtk_file.unstructured_grid.cells = {{0, 1, 3, 4}, {1, 2, 3}};
 //   vtk_file.unstructured_grid.cell_types = {VtkFile::CellType::VTK_QUAD,
 //                                            VtkFile::CellType::VTK_TRIANGLE};
 //   WriteToFile(vtk_file, "test.vtk");
 // }

}  // namespace lf::io::test
