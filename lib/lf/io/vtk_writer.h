/**
 * @file
 * @brief Declares the VtkWriter which can write vtk files that can be read
 *        by ParaView
 * @author Raffael Casagrande
 * @date   2018-07-14 07:13:49
 * @copyright MIT License
 */

#ifndef __3e48c7b32a034cb3be3dbca884ff4f6c
#define __3e48c7b32a034cb3be3dbca884ff4f6c

#include <Eigen/Eigen>
#include <string>
#include <vector>

namespace lf::io {

class VtkFile {
 public:
  using size_type = unsigned int;

  /// Nested types
  //////////////////////////////////////////////////////////////////////////////
  enum class Format { ASCII, BINARY };

  enum class CellType {
    VTK_VERTEX = 1,
    VTK_POLY_VERTEX = 2,
    VTK_LINE = 3,
    VTK_POLY_LINE = 4,
    VTK_TRIANGLE = 5,
    VTK_TRIANGLE_STRIP = 6,
    VTK_POLYGON = 7,
    VTK_PIXEL = 8,
    VTK_QUAD = 9,
    VTK_TETRA = 10,
    VTK_VOXEL = 11,
    VTK_HEXAHEDRON = 12,
    VTK_WEDGE = 13,
    VTK_PYRAMID = 14,
    VTK_QUADRATIC_EDGE = 21,
    VTK_QUADRATIC_TRIANGLE = 22,
    VTK_QUADRATIC_QUAD = 23,
    VTK_QUADRATIC_TETRA = 24,
    VTK_QUADRATIC_HEXAHEDRON = 25
  };

  class UnstructuredGrid {
   public:
    std::vector<Eigen::Vector3d> points;
    std::vector<std::vector<size_type>> cells;
    std::vector<CellType> cell_types;
  };

  /// Represents one set of attribute data (can be attached to points or cells)
  template <class T>
  class ScalarData {
   public:
    std::string data_name;
    std::vector<T> data;
    std::string lookup_table = "default";
  };

  class Attributes {
   public:
    std::vector<ScalarData<int>> scalar_int_data;
    std::vector<ScalarData<unsigned int>> scalar_unsigned_int_data;
    std::vector<ScalarData<long>> scalar_long_data;
    std::vector<ScalarData<unsigned long>> scalar_unsigned_long_data;
    std::vector<ScalarData<float>> scalar_float_data;
    std::vector<ScalarData<double>> scalar_double_data;
  };

  // Actual members
  //////////////////////////////////////////////////////////////////////////////

  /// The Vtk Header, at most 256 characters, no new lines characters!
  std::string header;
  /// The format of the file.
  Format format = Format::ASCII;

  /// Describes the nodes + cells
  UnstructuredGrid unstructured_grid;

  // Data that is attached to points
  Attributes point_data;
};

void WriteToFile(const VtkFile& vtk_file, const std::string& filename);

}  // namespace lf::io

#endif  // __3e48c7b32a034cb3be3dbca884ff4f6c
