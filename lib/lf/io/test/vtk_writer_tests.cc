/**
 * @file
 * @brief Test the vtk writer implementation.
 * @author Raffael Casagrande
 * @date   2018-07-14 07:52:25
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <complex>

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

TEST(lf_io_VtkWriter, vtkFilewriteOnlyMesh) {
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

TEST(lf_io_VtkWriter, twoElementMeshCodim0NoData) {
  auto reader = test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);

  // write mesh:
  VtkWriter writer(reader.mesh(), "two_element_no_data.vtk");
}

TEST(lf_io_VtkWriter, twoElementMeshCodim1NoData) {
  auto reader = test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);

  // write mesh:
  VtkWriter writer(reader.mesh(), "two_element_1d_nodata.vtk", 1);
}

template <class F>
struct MeshFunctionLambda {
  MeshFunctionLambda(F f) : f_(f) {}

  auto operator()(const mesh::Entity& e, const Eigen::MatrixXd& x) const {
    std::vector<std::invoke_result_t<F, const mesh::Entity&, Eigen::Matrix2d>>
        result(x.cols());
    for (int i = 0; i < x.cols(); ++i) {
      result[i] = f_(e, x.col(i));
    }
    return result;
  }

 private:
  F f_;
};

TEST(lf_io_VtkWriter, twoElementMeshCodim0AllData) {
  auto reader = test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);

  using lf::mesh::utils::make_LambdaMeshDataSet;

  // write mesh:
  VtkWriter writer(reader.mesh(), "two_element.vtk");
  writer.setBinary(true);

  writer.WriteGlobalData("global_zeros", std::vector<int>{0, 0, 0});
  writer.WriteGlobalData("global_ones", std::vector<double>{1., 1., 1.});
  writer.WriteGlobalData("global_twos", std::vector<float>{2., 2., 2.});

  auto zero = Eigen::VectorXd::Zero(0);
  auto zero2 = Eigen::Vector2d::Zero();

  writer.WritePointData(
      "uchar", *make_LambdaMeshDataSet([&](const auto& e) {
        return static_cast<unsigned char>(reader.mesh()->Index(e));
      }));
  writer.WritePointData(
      "uchar_mf",
      MeshFunctionLambda([&](const auto& e, const Eigen::MatrixXd& x) {
        return static_cast<unsigned char>(reader.mesh()->Index(e));
      }));
  writer.WritePointData(
      "char", *make_LambdaMeshDataSet(
                  [&](const auto& e) {
                    return static_cast<char>(reader.mesh()->Index(e));
                  },
                  [&](const auto& e) { return reader.mesh()->Index(e) < 3; }));
  writer.WritePointData(
      "char_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return static_cast<unsigned int>(reader.mesh()->Index(e));
      }));
  writer.WritePointData(
      "uint", *make_LambdaMeshDataSet([&](const auto& e) {
        return static_cast<unsigned int>(reader.mesh()->Index(e));
      }));
  writer.WritePointData(
      "uint_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return static_cast<unsigned int>(reader.mesh()->Index(e));
      }));
  writer.WritePointData("int", *make_LambdaMeshDataSet([&](const auto& e) {
                          return static_cast<int>(reader.mesh()->Index(e));
                        }));
  writer.WritePointData("int_mf",
                        MeshFunctionLambda([&](const auto& e, const auto& x) {
                          return static_cast<int>(reader.mesh()->Index(e));
                        }));
  writer.WritePointData("float", *make_LambdaMeshDataSet([&](const auto& e) {
                          return static_cast<float>(reader.mesh()->Index(e));
                        }));
  writer.WritePointData("float_mf",
                        MeshFunctionLambda([&](const auto& e, const auto& x) {
                          return static_cast<float>(reader.mesh()->Index(e));
                        }));
  writer.WritePointData("double", *make_LambdaMeshDataSet([&](const auto& e) {
                          return static_cast<double>(reader.mesh()->Index(e));
                        }));
  writer.WritePointData("double_mf",
                        MeshFunctionLambda([&](const auto& e, const auto& x) {
                          return static_cast<double>(reader.mesh()->Index(e));
                        }));
  writer.WritePointData("vector2d", *make_LambdaMeshDataSet([&](const auto& e) {
                          return Eigen::Vector2d(e.Geometry()->Global(zero));
                        }));
  writer.WritePointData("vector2d_mf",
                        MeshFunctionLambda([&](const auto& e, const auto& x) {
                          return Eigen::Vector2d(e.Geometry()->Global(zero));
                        }));
  writer.WritePointData(
      "vector2f", *make_LambdaMeshDataSet([&](const auto& e) {
        return Eigen::Vector2f(
            e.Geometry()->Global(zero).template cast<float>());
      }));
  writer.WritePointData(
      "vector2f_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return Eigen::Vector2f(
            e.Geometry()->Global(zero).template cast<float>());
      }));
  writer.WritePointData("vector3d", *make_LambdaMeshDataSet([&](const auto& e) {
                          return Eigen::Vector3d(0, 1, 2);
                        }));
  writer.WritePointData("vector3d_mf",
                        MeshFunctionLambda([&](const auto& e, const auto& x) {
                          return Eigen::Vector3d(0, 1, 2);
                        }));
  writer.WritePointData("vector3f", *make_LambdaMeshDataSet([&](const auto& e) {
                          return Eigen::Vector3f{0, 1, 2};
                        }));
  writer.WritePointData("vector3f_mf",
                        MeshFunctionLambda([&](const auto& e, const auto& x) {
                          return Eigen::Vector3f(0, 1, 2);
                        }));

  writer.WritePointData("vectorXd", *make_LambdaMeshDataSet(
                                        [&](const auto& e) -> Eigen::VectorXd {
                                          return Eigen::Vector2d(0, 1);
                                        }));
  writer.WritePointData(
      "vectorXd_mf",
      MeshFunctionLambda([&](const auto& e, const auto& x) -> Eigen::VectorXd {
        return Eigen::Vector2d(0, 1);
      }));
  writer.WritePointData("vectorXf", *make_LambdaMeshDataSet(
                                        [&](const auto& e) -> Eigen::VectorXf {
                                          return Eigen::Vector2f{0, 1};
                                        }));
  writer.WritePointData(
      "vectorXf_mf",
      MeshFunctionLambda([&](const auto& e, const auto& x) -> Eigen::VectorXf {
        return Eigen::Vector2f(0, 1);
      }));

  writer.WriteCellData(
      "uchar", *make_LambdaMeshDataSet([&](const auto& e) {
        return static_cast<unsigned char>(reader.mesh()->Index(e));
      }));
  writer.WriteCellData(
      "uchar_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return static_cast<unsigned char>(reader.mesh()->Index(e));
      }));
  writer.WriteCellData(
      "char", *make_LambdaMeshDataSet(
                  [&](const auto& e) {
                    return static_cast<char>(reader.mesh()->Index(e));
                  },
                  [&](const auto& e) { return reader.mesh()->Index(e) == 0; }));
  writer.WriteCellData(
      "uint", *make_LambdaMeshDataSet([&](const auto& e) {
        return static_cast<unsigned int>(reader.mesh()->Index(e));
      }));
  writer.WriteCellData(
      "uint_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return static_cast<unsigned int>(reader.mesh()->Index(e));
      }));
  writer.WriteCellData("int", *make_LambdaMeshDataSet([&](const auto& e) {
                         return static_cast<int>(reader.mesh()->Index(e));
                       }));
  writer.WriteCellData("int_mf",
                       MeshFunctionLambda([&](const auto& e, const auto& x) {
                         return static_cast<int>(reader.mesh()->Index(e));
                       }));
  writer.WriteCellData("float", *make_LambdaMeshDataSet([&](const auto& e) {
                         return static_cast<float>(reader.mesh()->Index(e));
                       }));
  writer.WriteCellData("float_mf",
                       MeshFunctionLambda([&](const auto& e, const auto& x) {
                         return static_cast<float>(reader.mesh()->Index(e));
                       }));
  writer.WriteCellData("double", *make_LambdaMeshDataSet([&](const auto& e) {
                         return static_cast<double>(reader.mesh()->Index(e));
                       }));
  writer.WriteCellData(
      "double_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return static_cast<unsigned int>(reader.mesh()->Index(e));
      }));
  writer.WriteCellData("vector2d", *make_LambdaMeshDataSet([&](const auto& e) {
                         return Eigen::Vector2d(e.Geometry()->Global(zero2));
                       }));
  writer.WriteCellData("vector2d_mf",
                       MeshFunctionLambda([&](const auto& e, const auto& x) {
                         return Eigen::Vector2d(e.Geometry()->Global(zero2));
                       }));
  writer.WriteCellData(
      "vector2f", *make_LambdaMeshDataSet([&](const auto& e) {
        return Eigen::Vector2f(
            e.Geometry()->Global(zero2).template cast<float>());
      }));
  writer.WriteCellData(
      "vector2f_mf", MeshFunctionLambda([&](const auto& e, const auto& x) {
        return Eigen::Vector2f(
            e.Geometry()->Global(zero2).template cast<float>());
      }));
  writer.WriteCellData("vector3d", *make_LambdaMeshDataSet([&](const auto& e) {
                         return Eigen::Vector3d(0, 1, 2);
                       }));
  writer.WriteCellData("vector3d_mf",
                       MeshFunctionLambda([&](const auto& e, const auto& x) {
                         return Eigen::Vector3d(0, 1, 2);
                       }));
  writer.WriteCellData("vector3f", *make_LambdaMeshDataSet([&](const auto& e) {
                         return Eigen::Vector3f{0, 1, 2};
                       }));
  writer.WriteCellData("vector3f_mf",
                       MeshFunctionLambda([&](const auto& e, const auto& x) {
                         return Eigen::Vector3f{0, 1, 2};
                       }));
  writer.WriteCellData("vectorXd", *make_LambdaMeshDataSet(
                                       [&](const auto& e) -> Eigen::VectorXd {
                                         return Eigen::Vector3d(0, 1, 2);
                                       }));
  writer.WriteCellData(
      "vectorXd_mf",
      MeshFunctionLambda([&](const auto& e, const auto& x) -> Eigen::VectorXd {
        return Eigen::Vector3d(0, 1, 2);
      }));

  writer.WriteCellData("vectorXf", *make_LambdaMeshDataSet(
                                       [&](const auto& e) -> Eigen::VectorXf {
                                         return Eigen::Vector3f(0, 1, 2);
                                       }));
  writer.WriteCellData(
      "vectorXf_mf",
      MeshFunctionLambda([&](const auto& e, const auto& x) -> Eigen::VectorXf {
        return Eigen::Vector3f(0, 1, 2);
      }));

  // try to write data with a name that is already used:
  EXPECT_THROW(
      writer.WritePointData(
          "char", *make_LambdaMeshDataSet([](const auto& e) { return 1; })),
      base::LfException);
  EXPECT_THROW(
      writer.WriteCellData(
          "char", *make_LambdaMeshDataSet([](const auto& e) { return 1; })),
      base::LfException);
  EXPECT_THROW(writer.WriteGlobalData("global_zeros", std::vector<int>{0}),
               base::LfException);

  // try to write data with a name that contains spaces:
  EXPECT_THROW(
      writer.WritePointData(
          "h w", *make_LambdaMeshDataSet([](const auto& e) { return 1; })),
      base::LfException);
  EXPECT_THROW(
      writer.WriteCellData(
          "h w", *make_LambdaMeshDataSet([](const auto& e) { return 1; })),
      base::LfException);
  EXPECT_THROW(writer.WriteGlobalData("h w", std::vector<int>{0}),
               base::LfException);
}

void writeSomeData(VtkWriter& vtk) {
  auto scalarfn = [&](const Eigen::Vector2d& x) {
    return std::pow(std::complex<double>(x[0], x[1]), 6).real();
  };

  vtk.WritePointData("linear",
                     mesh::utils::MeshFunctionGlobal(
                         [](const Eigen::Vector2d& x) { return x[0] + x[1]; }));

  vtk.WritePointData("unsigned_char",
                     mesh::utils::MeshFunctionGlobal(
                         [&](const Eigen::Vector2d& x) -> unsigned char {
                           return scalarfn(x);
                         }));
  vtk.WritePointData("char", mesh::utils::MeshFunctionGlobal(
                                 [&](const Eigen::Vector2d& x) -> char {
                                   return scalarfn(x);
                                 }));
  vtk.WritePointData("unsigned", mesh::utils::MeshFunctionGlobal(
                                     [&](const Eigen::Vector2d& x) -> unsigned {
                                       return scalarfn(x);
                                     }));
  vtk.WritePointData(
      "int", mesh::utils::MeshFunctionGlobal(
                 [&](const Eigen::Vector2d& x) -> int { return scalarfn(x); }));
  vtk.WritePointData("float", mesh::utils::MeshFunctionGlobal(
                                  [&](const Eigen::Vector2d& x) -> float {
                                    return scalarfn(x);
                                  }));
  vtk.WritePointData("double", mesh::utils::MeshFunctionGlobal(scalarfn));
  vtk.WritePointData("vec2_float",
                     mesh::utils::MeshFunctionGlobal(
                         [](const Eigen::Vector2d& x) -> Eigen::Vector2f {
                           return x.cast<float>();
                         }));
  vtk.WritePointData("vec2_double",
                     mesh::utils::MeshFunctionGlobal(
                         [](const Eigen::Vector2d& x) { return x; }));
  vtk.WritePointData(
      "vec3_float",
      mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return Eigen::Vector3f(x[0], x[1], std::sin(x[0]));
      }));
  vtk.WritePointData(
      "vec3_double",
      mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return Eigen::Vector3d(x[0], x[1], std::sin(x[0]));
      }));
  vtk.WritePointData(
      "vecd_float",
      mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return Eigen::VectorXf(x.cast<float>());
      }));
  vtk.WritePointData("vecd_double", mesh::utils::MeshFunctionGlobal(
                                        [](const Eigen::Vector2d& x) {
                                          return Eigen::VectorXd(x);
                                        }));
}

TEST(lf_io_VtkWriter, circle2ndOrderQuad) {
  auto reader = test_utils::getGmshReader("circle_second_order_quad.msh", 2);
  VtkWriter vtk(reader.mesh(), "circle_second_order_quad.vtk", 0, 2);
  writeSomeData(vtk);
  VtkWriter vtk2(reader.mesh(), "circle_second_order_quad_codim1.vtk", 1, 5);
  writeSomeData(vtk2);

  auto readerTria = test_utils::getGmshReader("circle_second_order.msh", 2);
  VtkWriter vtkt(readerTria.mesh(), "circle_second_order_tria.vtk", 0, 2);
  writeSomeData(vtkt);
}
}  // namespace lf::io::test
