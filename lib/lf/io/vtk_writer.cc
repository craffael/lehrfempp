/**
 * @file
 * @brief Implementation of VtkWriter
 * @author Raffael Casagrande
 * @date   2018-07-14 07:18:02
 * @copyright MIT License
 */

#include "vtk_writer.h"
#include <lf/base/base.h>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/phoenix/scope/let.hpp>
#include <boost/spirit/include/karma.hpp>
#include <fstream>
#include "eigen_fusion_adapter.h"

template <class RESULT_TYPE, class LAMBDA>
class VariantVisitor {
 public:
  using result_type = RESULT_TYPE;

  template <class T>
  result_type operator()(T&& arg) {
    return lambda_(std::forward<T>(arg));
  }

  explicit VariantVisitor(LAMBDA lambda) : lambda_(lambda) {}

 private:
  LAMBDA lambda_;
};

template <class RESULT_TYPE, class LAMBDA>
VariantVisitor<RESULT_TYPE, LAMBDA> make_VariantVisitor(LAMBDA l) {
  return VariantVisitor<RESULT_TYPE, LAMBDA>(l);
}

// clang-format off
BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile::UnstructuredGrid,
    (std::vector<Eigen::Vector3f>, points)
    (std::vector<std::vector<unsigned int>>, cells)
    (std::vector<lf::io::VtkFile::CellType>, cell_types)
);

BOOST_FUSION_ADAPT_TPL_STRUCT((T), (lf::io::VtkFile::FieldDataArray)(T),
                              (std::string, name),
                              (auto, data)
);

BOOST_FUSION_ADAPT_TPL_STRUCT((T), (lf::io::VtkFile::ScalarData)(T),
    (std::string, name)
    (std::string, lookup_table)
    (auto, data)
);

BOOST_FUSION_ADAPT_TPL_STRUCT((T), (lf::io::VtkFile::VectorData)(T),
    (std::string, name)
    (auto, data)
);

// BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile::Attributes,
//                           (auto, data)
//                           );

 // BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile::Attributes,
 //     (std::vector<lf::io::VtkFile::ScalarData<int>>, scalar_int_data)
 //     (std::vector<lf::io::VtkFile::ScalarData<unsigned int>>, scalar_unsigned_int_data)
 //     (std::vector<lf::io::VtkFile::ScalarData<long>>, scalar_long_data)
 //     (std::vector<lf::io::VtkFile::ScalarData<unsigned long>>, scalar_unsigned_long_data)
 //     (std::vector<lf::io::VtkFile::ScalarData<float>>, scalar_float_data)
 //     (std::vector<lf::io::VtkFile::ScalarData<double>>, scalar_double_data)
 // );

BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile,
    (std::string, header)
    (lf::io::VtkFile::Format, format)
    (lf::io::VtkFile::UnstructuredGrid, unstructured_grid)
    (lf::io::VtkFile::FieldData, field_data)
    (lf::io::VtkFile::Attributes, point_data)
    (lf::io::VtkFile::Attributes, cell_data)
);

// clang-format on

namespace /*anonymous*/ {
void cell_list_size(unsigned int& result,  // NOLINT
                    const std::vector<std::vector<unsigned int>>& cells) {
  result = cells.size();
  for (auto& v : cells) {
    result += v.size();
  }
}
}  // namespace

BOOST_PHOENIX_ADAPT_FUNCTION(void, CellListSize, cell_list_size, 2);

namespace lf::io {
std::ostream& operator<<(std::ostream& stream, const VtkFile::Format& f) {
  switch (f) {
    case VtkFile::Format::ASCII:
      stream << "ASCII";
      break;
    case VtkFile::Format::BINARY:
      stream << "BINARY";
      break;
    default:
      throw base::LfException("Unknown VtkFileFormat.");
  }
  return stream;
}

namespace /*anonymous */ {
namespace karma = boost::spirit::karma;
template <class Iterator, bool BINARY>
struct VtkGrammar : karma::grammar<Iterator, VtkFile()> {
  // template <class T>
  // karma::rule<Iterator, std::vector<T>()> vector_size() {
  //   return karma::int_[karma::_1 = boost::phoenix::size(karma::_val)];
  // }

  VtkGrammar() : VtkGrammar::base_type(root) {
    using boost::phoenix::at_c;
    using boost::phoenix::size;
    using karma::_1;
    using karma::_a;
    using karma::_r1;
    using karma::_val;
    using karma::eol;
    using karma::eps;
    using karma::lit;

    if (BINARY) {
      points_vector %= *(karma::big_bin_float << karma::big_bin_float
                                              << karma::big_bin_float)
                       << eol;
    } else {
      points_vector %= *(karma::float_ << ' ' << karma::float_ << ' '
                                       << karma::float_ << karma::eol);
    }

    points = lit("POINTS ")
             << karma::uint_[karma::_1 = boost::phoenix::size(karma::_val)]
             << lit(" float") << karma::eol
             << points_vector[karma::_1 = karma::_val];

    // cells_vector %= (karma::uint_ % " ") << karma::eol;
    if (BINARY) {
      cells_vector %= karma::big_dword(boost::phoenix::size(karma::_val))
                      << *(karma::big_dword);
    } else {
      // cells_vector %= karma::duplicate
      //     [karma::uint_[karma::_1 = boost::phoenix::size(karma::_val)]
      //      << lit(' ') << (karma::uint_ % " ") << karma::eol];
      cells_vector %= lit(size(_val))
                      << lit(' ') << (karma::uint_ % " ") << karma::eol;
    }
    cells = lit("CELLS ")
            << karma::uint_[karma::_1 = boost::phoenix::size(karma::_val)]
            << lit(' ') << karma::uint_[CellListSize(karma::_1, karma::_val)]
            << karma::eol << (*cells_vector)[karma::_1 = karma::_val] << eol;

    if (BINARY) {
      cell_type =
          karma::big_dword[karma::_1 =
                               boost::phoenix::static_cast_<int>(karma::_val)];
    } else {
      cell_type = karma::int_[karma::_1 = boost::phoenix::static_cast_<int>(
                                  karma::_val)]
                  << karma::eol;
    }
    cell_types = lit("CELL_TYPES ")
                 << karma::uint_[karma::_1 = boost::phoenix::size(karma::_val)]
                 << karma::eol << (*cell_type)[karma::_1 = karma::_val];

    unstructured_grid %= lit("DATASET UNSTRUCTURED_GRID")
                         << karma::eol << points << cells << cell_types;

    // Field data
    if (BINARY) {
      field_data_int %= karma::string << lit(" 1 ") << lit(size(at_c<1>(_val)))
                                      << lit(" int") << eol << *karma::big_dword
                                      << eol;
      field_data_float %= karma::string
                          << lit(" 1 ") << lit(size(at_c<1>(_val)))
                          << lit(" float") << eol << *(karma::big_bin_float)
                          << eol;

      field_data_double %= karma::string
                           << lit(" 1 ") << lit(size(at_c<1>(_val)))
                           << lit(" double") << eol << *karma::big_bin_double
                           << eol;
    } else {
      field_data_int %= karma::string << lit(" 1 ") << lit(size(at_c<1>(_val)))
                                      << lit(" int") << eol
                                      << (karma::int_ % lit(' ')) << eol;
      field_data_float %= karma::string
                          << lit(" 1 ") << lit(size(at_c<1>(_val)))
                          << lit(" float") << eol << (karma::float_ % lit(' '))
                          << eol;

      field_data_double %= karma::string << lit(" 1 ")
                                         << lit(size(at_c<1>(_val)))
                                         << lit(" double") << eol
                                         << (karma::double_ % lit(' ')) << eol;
    }
    field_data %= lit("FIELD FieldData ")
                  << lit(size(_val)) << eol
                  << *(field_data_int | field_data_float | field_data_double);

    // Point/Cell data
    if (BINARY) {
      scalar_data_char %=
          (lit("SCALARS ") << karma::string << lit(" char 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::byte_)
          << eol;

      scalar_data_uchar %=
          (lit("SCALARS ") << karma::string << lit(" char 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::byte_)
          << eol;

      scalar_data_short %=
          (lit("SCALARS ") << karma::string << lit(" short 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::big_word)
          << eol;

      scalar_data_ushort %=
          (lit("SCALARS ") << karma::string << lit(" unsigned_short 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::big_word)
          << eol;

      scalar_data_int %=
          (lit("SCALARS ") << karma::string << lit(" int 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::big_dword)
          << eol;
      scalar_data_unsigned_int %=
          (lit("SCALARS ") << karma::string << lit(" unsigned_int 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::big_dword)
          << eol;

      // NOLINTNEXTLINE
      if constexpr (sizeof(long) == 8) {
        scalar_data_long %=
            (lit("SCALARS ")
             << karma::string << lit(" long 1") << eol << lit("LOOKUP_TABLE ")
             << karma::string << eol << *karma::big_qword)
            << eol;
        scalar_data_unsigned_long %=
            (lit("SCALARS ") << karma::string << lit(" unsigned_long 1") << eol
                             << lit("LOOKUP_TABLE ") << karma::string << eol
                             << *karma::big_qword)
            << eol;
      } else {
        scalar_data_long %=
            (lit("SCALARS ")
             << karma::string << lit(" long 1") << eol << lit("LOOKUP_TABLE ")
             << karma::string << eol << *karma::big_dword)
            << eol;
        scalar_data_unsigned_long %=
            (lit("SCALARS ") << karma::string << lit(" unsigned_long 1") << eol
                             << lit("LOOKUP_TABLE ") << karma::string << eol
                             << *karma::big_dword)
            << eol;
      }

      scalar_data_float %=
          (lit("SCALARS ") << karma::string << lit(" float 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::big_bin_float)
          << eol;

      scalar_data_double %=
          (lit("SCALARS ") << karma::string << lit(" double 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << *karma::big_bin_double)
          << eol;

      vector_data_float %= lit("VECTORS ")
                           << karma::string << lit(" float") << eol
                           << *(karma::big_bin_float << karma::big_bin_float
                                                     << karma::big_bin_float)
                           << eol;

      vector_data_double %= lit("VECTORS ")
                            << karma::string << lit(" double") << eol
                            << *(karma::big_bin_double << karma::big_bin_double
                                                       << karma::big_bin_double)
                            << eol;
    } else {
      scalar_data_char %=
          (lit("SCALARS ") << karma::string << lit(" char 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::short_ % eol)))
          << eol;

      scalar_data_uchar %=
          (lit("SCALARS ") << karma::string << lit(" char 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::ushort_ % eol)))
          << eol;

      scalar_data_short %=
          (lit("SCALARS ") << karma::string << lit(" short 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::short_ % eol)))
          << eol;

      scalar_data_ushort %=
          (lit("SCALARS ") << karma::string << lit(" unsigned_short 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::ushort_ % eol)))
          << eol;

      scalar_data_int %=
          (lit("SCALARS ") << karma::string << lit(" int 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::int_ % eol)))
          << eol;
      scalar_data_unsigned_int %=
          (lit("SCALARS ") << karma::string << lit(" unsigned_int 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::uint_ % eol)))
          << eol;

      scalar_data_long %=
          (lit("SCALARS ") << karma::string << lit(" long 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::long_ % eol)))
          << eol;

      scalar_data_unsigned_long %=
          (lit("SCALARS ") << karma::string << lit(" unsigned_long 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::ulong_ % eol)))
          << eol;

      scalar_data_float %=
          (lit("SCALARS ") << karma::string << lit(" float 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::float_ % eol)))
          << eol;

      scalar_data_double %=
          (lit("SCALARS ") << karma::string << lit(" double 1") << eol
                           << lit("LOOKUP_TABLE ") << karma::string << eol
                           << (*(karma::double_ % eol)))
          << eol;

      vector_data_float %= lit("VECTORS ")
                           << karma::string << lit(" float") << eol
                           << *(karma::float_ << lit(' ') << karma::float_
                                              << lit(' ') << karma::float_
                                              << eol);

      vector_data_double %= lit("VECTORS ")
                            << karma::string << lit(" double") << eol
                            << *(karma::double_ << lit(' ') << karma::double_
                                                << lit(' ') << karma::double_
                                                << eol);
    }

    attributes %=
        (*(scalar_data_char | scalar_data_uchar | scalar_data_short |
           scalar_data_ushort | scalar_data_int | scalar_data_unsigned_int |
           scalar_data_long | scalar_data_unsigned_long | scalar_data_float |
           scalar_data_double | vector_data_float | vector_data_double));

    root %= lit("# vtk DataFile Version 3.0")
            << karma::eol << karma::string << karma::eol << karma::stream
            << karma::eol << unstructured_grid << karma::eol << field_data
            << lit("POINT_DATA ") << lit(size(at_c<0>(at_c<2>(_val)))) << eol
            << attributes << lit("CELL_DATA ")
            << lit(size(at_c<1>(at_c<2>(_val)))) << eol << attributes;
  }

  karma::rule<Iterator, VtkFile()> root;
  karma::rule<Iterator, VtkFile()> header;
  karma::rule<Iterator, VtkFile::UnstructuredGrid()> unstructured_grid;
  karma::rule<Iterator, std::vector<Eigen::Vector3f>()> points;
  karma::rule<Iterator, std::vector<Eigen::Vector3f>()> points_vector;
  karma::rule<Iterator, std::vector<std::vector<unsigned int>>> cells;
  karma::rule<Iterator, std::vector<unsigned int>> cells_vector;
  karma::rule<Iterator, std::vector<VtkFile::CellType>> cell_types;
  karma::rule<Iterator, VtkFile::CellType> cell_type;

  karma::rule<Iterator, VtkFile::FieldData()> field_data;
  karma::rule<Iterator, VtkFile::FieldDataArray<int>()> field_data_int;
  karma::rule<Iterator, VtkFile::FieldDataArray<float>()> field_data_float;
  karma::rule<Iterator, VtkFile::FieldDataArray<double>()> field_data_double;

  karma::rule<Iterator, VtkFile::Attributes()> attributes;

  karma::rule<Iterator, VtkFile::ScalarData<char>()> scalar_data_char;
  karma::rule<Iterator, VtkFile::ScalarData<unsigned char>()> scalar_data_uchar;
  karma::rule<Iterator, VtkFile::ScalarData<short>()> scalar_data_short;
  karma::rule<Iterator, VtkFile::ScalarData<unsigned short>()>
      scalar_data_ushort;
  karma::rule<Iterator, VtkFile::ScalarData<int>()> scalar_data_int;
  karma::rule<Iterator, VtkFile::ScalarData<unsigned int>()>
      scalar_data_unsigned_int;
  karma::rule<Iterator, VtkFile::ScalarData<long>()> scalar_data_long;
  karma::rule<Iterator, VtkFile::ScalarData<unsigned long>()>
      scalar_data_unsigned_long;
  karma::rule<Iterator, VtkFile::ScalarData<float>()> scalar_data_float;
  karma::rule<Iterator, VtkFile::ScalarData<double>()> scalar_data_double;
  karma::rule<Iterator, VtkFile::VectorData<float>()> vector_data_float;
  karma::rule<Iterator, VtkFile::VectorData<double>()> vector_data_double;
};

/// Perform some checks to make sure the vtk file is valid.
void ValidateVtkFile(const VtkFile& vtk_file) {
  // check that header doesn't contain a new line character:
  if (vtk_file.header.find('\n') != std::string::npos) {
    throw base::LfException(
        "Header of vtk file contains a new line character, this is not "
        "allowed");
  }
  if (vtk_file.header.size() > 256) {
    throw base::LfException("Header of vtk file is longer than 256 characters");
  }
  if (vtk_file.unstructured_grid.cell_types.size() !=
      vtk_file.unstructured_grid.cells.size()) {
    throw base::LfException(
        "Mismatch of size of cell_types and cells in VtkFile.");
  }
  for (auto& d : vtk_file.point_data) {
    boost::apply_visitor(
        [&](auto e) {
          if (e.data.size() != vtk_file.unstructured_grid.points.size()) {
            throw base::LfException("Length of dataset " + e.name +
                                    " does not match the number of points.");
          }
        },
        d);
  }
  for (auto& d : vtk_file.cell_data) {
    boost::apply_visitor(
        [&](auto e) {
          if (e.data.size() != vtk_file.unstructured_grid.cells.size()) {
            throw base::LfException("Length of dataset " + e.name +
                                    " does not match the number of cells.");
          }
        },
        d);
  }
}
}  // namespace

void WriteToFile(const VtkFile& vtk_file, const std::string& filename) {
  std::ofstream file(filename, std::ios_base::out | std::ios_base::binary |
                                   std::ios_base::trunc);
  ValidateVtkFile(vtk_file);

  if (!file.is_open()) {
    throw base::LfException("Could not open file " + filename +
                            " for writing.");
  }
  karma::ostream_iterator<char> outit(file);

  bool result = false;
  if (vtk_file.format == VtkFile::Format::BINARY) {
    VtkGrammar<decltype(outit), true> grammar{};
    result = karma::generate(outit, grammar, vtk_file);
  } else {
    VtkGrammar<decltype(outit), false> grammar{};
    result = karma::generate(outit, grammar, vtk_file);
  }

  file.close();
  if (!result) {
    throw base::LfException("Karma error");
  }
}

VtkWriter::VtkWriter(std::shared_ptr<mesh::Mesh> mesh, std::string filename,
                     dim_t codim)
    : mesh_(std::move(mesh)), filename_(filename), codim_(codim) {
  auto dim_mesh = mesh_->DimMesh();
  auto dim_world = mesh_->DimWorld();
  LF_ASSERT_MSG(dim_world > 0 && dim_world <= 4,
                "VtkWriter supports only dim_world = 1,2 or 3");
  LF_ASSERT_MSG(codim >= 0 && codim < dim_mesh, "codim out of bounds.");

  // insert nodes:
  vtk_file_.unstructured_grid.points.resize(mesh_->Size(dim_mesh));
  Eigen::Matrix<double, 0, 1> zero;
  for (auto& p : mesh_->Entities(dim_mesh)) {
    auto index = mesh_->Index(p);
    Eigen::Vector3f coord;
    if (dim_world == 1) {
      coord(0) = p.Geometry()->Global(zero)(0);
      coord(1) = 0.f;
      coord(2) = 0.f;
    } else if (dim_world == 2) {
      coord.topRows<2>() = p.Geometry()->Global(zero).cast<float>();
      coord(2) = 0.f;
    } else {
      coord = p.Geometry()->Global(zero).cast<float>();
    }

    vtk_file_.unstructured_grid.points[index] = std::move(coord);
  }

  // insert elements:
  vtk_file_.unstructured_grid.cells.resize(mesh_->Size(codim));
  vtk_file_.unstructured_grid.cell_types.resize(mesh_->Size(codim));
  for (auto& e : mesh_->Entities(codim)) {
    auto index = mesh_->Index(e);
    auto ref_el = e.RefEl();
    auto& node_indices = vtk_file_.unstructured_grid.cells[index];
    node_indices.reserve(ref_el.NumNodes());
    for (auto& p : e.SubEntities(dim_mesh - codim)) {
      node_indices.push_back(mesh_->Index(p));
    }

    if (ref_el == base::RefEl::kSegment()) {
      vtk_file_.unstructured_grid.cell_types[index] =
          VtkFile::CellType::VTK_LINE;
    } else if (ref_el == base::RefEl::kTria()) {
      vtk_file_.unstructured_grid.cell_types[index] =
          VtkFile::CellType::VTK_TRIANGLE;
    } else if (ref_el == base::RefEl::kQuad()) {
      vtk_file_.unstructured_grid.cell_types[index] =
          VtkFile::CellType::VTK_QUAD;
    }
  }
}

void VtkWriter::WritePointData(
    const std::string& name, const mesh::utils::MeshDataSet<unsigned char>& mds,
    unsigned char undefined_value) {
  WriteScalarPointData(name, mds, undefined_value);
}

void VtkWriter::WritePointData(const std::string& name,
                               const mesh::utils::MeshDataSet<char>& mds,
                               char undefined_value) {
  WriteScalarPointData(name, mds, undefined_value);
}

void VtkWriter::WritePointData(const std::string& name,
                               const mesh::utils::MeshDataSet<unsigned>& mds,
                               unsigned undefined_value) {
  WriteScalarPointData(name, mds, undefined_value);
}

void VtkWriter::WritePointData(const std::string& name,
                               const mesh::utils::MeshDataSet<int>& mds,
                               int undefined_value) {
  WriteScalarPointData(name, mds, undefined_value);
}

void VtkWriter::WritePointData(const std::string& name,
                               const mesh::utils::MeshDataSet<float>& mds,
                               float undefined_value) {
  WriteScalarPointData(name, mds, undefined_value);
}

void VtkWriter::WritePointData(const std::string& name,
                               const mesh::utils::MeshDataSet<double>& mds,
                               double undefined_value) {
  WriteScalarPointData(name, mds, undefined_value);
}

void VtkWriter::WritePointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector2d>& mds,
    const Eigen::Vector2d& undefined_value) {
  WriteVectorPointData<2, double>(name, mds, undefined_value);
}

void VtkWriter::WritePointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector2f>& mds,
    const Eigen::Vector2f& undefined_value) {
  WriteVectorPointData<2, float>(name, mds, undefined_value);
}

void VtkWriter::WritePointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector3d>& mds,
    const Eigen::Vector3d& undefined_value) {
  WriteVectorPointData<3, double>(name, mds, undefined_value);
}

void VtkWriter::WritePointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector3f>& mds,
    const Eigen::Vector3f& undefined_value) {
  WriteVectorPointData<3, float>(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(
    const std::string& name, const mesh::utils::MeshDataSet<unsigned char>& mds,
    unsigned char undefined_value) {
  WriteScalarCellData(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(const std::string& name,
                              const mesh::utils::MeshDataSet<char>& mds,
                              char undefined_value) {
  WriteScalarCellData(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(const std::string& name,
                              const mesh::utils::MeshDataSet<unsigned>& mds,
                              unsigned undefined_value) {
  WriteScalarCellData(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(const std::string& name,
                              const mesh::utils::MeshDataSet<int>& mds,
                              int undefined_value) {
  WriteScalarCellData(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(const std::string& name,
                              const mesh::utils::MeshDataSet<float>& mds,
                              float undefined_value) {
  WriteScalarCellData(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(const std::string& name,
                              const mesh::utils::MeshDataSet<double>& mds,
                              double undefined_value) {
  WriteScalarCellData(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector2d>& mds,
    const Eigen::Vector2d& undefined_value) {
  WriteVectorCellData<2, double>(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector2f>& mds,
    const Eigen::Vector2f& undefined_value) {
  WriteVectorCellData<2, float>(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector3d>& mds,
    const Eigen::Vector3d& undefined_value) {
  WriteVectorCellData<3, double>(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Vector3f>& mds,
    const Eigen::Vector3f& undefined_value) {
  WriteVectorCellData<3, float>(name, mds, undefined_value);
}

void VtkWriter::WriteGlobalData(const std::string& name,
                                std::vector<int> data) {
  WriteFieldData(name, std::move(data));
}

void VtkWriter::WriteGlobalData(const std::string& name,
                                std::vector<float> data) {
  WriteFieldData(name, std::move(data));
}

void VtkWriter::WriteGlobalData(const std::string& name,
                                std::vector<double> data) {
  WriteFieldData(name, std::move(data));
}

template <class DATA>
void CheckAttributeSetName(const DATA& data, const std::string& name) {
  if (std::find_if(data.begin(), data.end(), [&](auto& d) {
        return boost::apply_visitor([&](auto&& d2) { return d2.name; }, d) ==
               name;
      }) != data.end()) {
    throw base::LfException(
        "There is already another Point/Cell Attribute Set with the name " +
        name);
  }
  if (name.find(' ') != std::string::npos) {
    throw base::LfException(
        "The name of the attribute set cannot contain spaces!");
  }
}

template <class T>
void VtkWriter::WriteScalarPointData(const std::string& name,
                                     const mesh::utils::MeshDataSet<T>& mds,
                                     T undefined_value) {
  CheckAttributeSetName(vtk_file_.point_data, name);
  VtkFile::ScalarData<T> data{};
  data.data.resize(mesh_->Size(mesh_->DimMesh()));
  data.name = name;
  for (auto& p : mesh_->Entities(mesh_->DimMesh())) {
    if (mds.DefinedOn(p)) {
      data.data[mesh_->Index(p)] = mds(p);
    } else {
      data.data[mesh_->Index(p)] = undefined_value;
    }
  }
  vtk_file_.point_data.push_back(std::move(data));
}

template <int ROWS, class T>
void PadWithZeros(Eigen::Matrix<T, 3, 1>& out,
                  const Eigen::Matrix<T, ROWS, 1>& in) {
  if constexpr (ROWS == 2) {  // NOLINT
    out.template block<2, 1>(0, 0) = in;
    out(2) = T(0);
  } else if constexpr (ROWS == 3) {  // NOLINT
    out = in;
  } else if constexpr (ROWS == Eigen::Dynamic) {  // NOLINT
    out(2) = T(0);
    out.topRows(in.rows()) = in;
  }
}

template <int ROWS, class T>
void VtkWriter::WriteVectorPointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Matrix<T, ROWS, 1>>& mds,
    const Eigen::Matrix<T, ROWS, 1>& undefined_value) {
  CheckAttributeSetName(vtk_file_.point_data, name);
  VtkFile::VectorData<T> data{};
  data.data.resize(mesh_->Size(mesh_->DimMesh()));
  data.name = name;
  Eigen::Matrix<T, 3, 1> undefined_value_padded;
  PadWithZeros<ROWS, T>(undefined_value_padded, undefined_value);

  for (auto& p : mesh_->Entities(mesh_->DimMesh())) {
    if (mds.DefinedOn(p)) {
      PadWithZeros<ROWS, T>(data.data[mesh_->Index(p)], mds(p));
    } else {
      data.data[mesh_->Index(p)] = undefined_value_padded;
    }
  }
  vtk_file_.point_data.push_back(std::move(data));
}

template <class T>
void VtkWriter::WriteScalarCellData(const std::string& name,
                                    const mesh::utils::MeshDataSet<T>& mds,
                                    T undefined_value) {
  CheckAttributeSetName(vtk_file_.cell_data, name);
  VtkFile::ScalarData<T> data{};
  data.data.resize(mesh_->Size(codim_));
  data.name = name;
  for (auto& e : mesh_->Entities(codim_)) {
    if (mds.DefinedOn(e)) {
      data.data[mesh_->Index(e)] = mds(e);
    } else {
      data.data[mesh_->Index(e)] = undefined_value;
    }
  }
  vtk_file_.cell_data.push_back(std::move(data));
}

template <int ROWS, class T>
void VtkWriter::WriteVectorCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Matrix<T, ROWS, 1>>& mds,
    const Eigen::Matrix<T, ROWS, 1>& undefined_value) {
  CheckAttributeSetName(vtk_file_.cell_data, name);
  VtkFile::VectorData<T> data{};
  data.data.resize(mesh_->Size(codim_));
  data.name = name;
  Eigen::Matrix<T, 3, 1> undefined_value_padded;
  PadWithZeros<ROWS, T>(undefined_value_padded, undefined_value);

  for (auto& p : mesh_->Entities(codim_)) {
    if (mds.DefinedOn(p)) {
      PadWithZeros<ROWS, T>(data.data[mesh_->Index(p)], mds(p));
    } else {
      data.data[mesh_->Index(p)] = undefined_value_padded;
    }
  }
  vtk_file_.cell_data.push_back(std::move(data));
}

template <class T>
void VtkWriter::WriteFieldData(const std::string& name, std::vector<T> data) {
  CheckAttributeSetName(vtk_file_.field_data, name);
  VtkFile::FieldDataArray<T> vtk_data{};
  vtk_data.name = name;
  vtk_data.data = std::move(data);
  vtk_file_.field_data.push_back(std::move(vtk_data));
}

}  // namespace lf::io
