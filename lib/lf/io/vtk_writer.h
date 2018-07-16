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
#include <boost/variant/variant.hpp>
#include <string>
#include <variant>
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

  template <class T>
  class FieldDataArray {
   public:
    std::string name;
    std::vector<T> data;

    FieldDataArray() {}

    explicit FieldDataArray(std::string name, std::vector<T> data)
        : name(name), data(data) {}
  };

  /// Represents one set of attribute data (can be attached to points or cells)
  template <class T>
  class ScalarData {
   public:
    std::string data_name;
    std::vector<T> data;
    std::string lookup_table = "default";

    ScalarData() {}

    explicit ScalarData(std::string data_name, std::vector<T> data,
                        std::string lookup_table = "default")
        : data_name(data_name),
          data(std::move(data)),
          lookup_table(lookup_table) {}
  };

  template <class T>
  class VectorData {
   public:
    std::string name;
    std::vector<Eigen::Matrix<T, 3, 1>> data;

    VectorData() {}

    explicit VectorData(std::string name,
                        std::vector<Eigen::Matrix<T, 3, 1>> data)
        : name(name), data(data) {}
  };

  using Attributes = std::vector<boost::variant<
      ScalarData<char>, ScalarData<unsigned char>, ScalarData<short>,
      ScalarData<unsigned short>, ScalarData<int>, ScalarData<unsigned int>,
      ScalarData<long>, ScalarData<unsigned long>, ScalarData<float>,
      ScalarData<double>, VectorData<float>, VectorData<double>>>;

  using FieldData =
      std::vector<boost::variant<FieldDataArray<int>, FieldDataArray<float>,
                                 FieldDataArray<double>>>;

  // Actual members
  //////////////////////////////////////////////////////////////////////////////

  /// The Vtk Header, at most 256 characters, no new lines characters!
  std::string header;
  /// The format of the file.
  Format format = Format::ASCII;

  /// Describes the nodes + cells
  UnstructuredGrid unstructured_grid;

  FieldData field_data;

  // Data that is attached to points
  Attributes point_data;

  // Data that is attached to the cells of the mesh
  Attributes cell_data;
};

void WriteToFile(const VtkFile& vtk_file, const std::string& filename);

}  // namespace lf::io

#endif  // __3e48c7b32a034cb3be3dbca884ff4f6c
