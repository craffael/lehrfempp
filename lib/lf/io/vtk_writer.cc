/**
 * @file
 * @brief Implementation of VtkWriter
 * @author Raffael Casagrande
 * @date   2018-07-14 07:18:02
 * @copyright MIT License
 */

#include "vtk_writer.h"

#include <lf/base/base.h>

#include <Eigen/Eigen>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/phoenix/phoenix.hpp>
#include <boost/phoenix/scope/let.hpp>
#include <boost/spirit/include/karma.hpp>
#include <fstream>
#include <unsupported/Eigen/KroneckerProduct>

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
  for (const auto& v : cells) {
    result += v.size();
  }
}
}  // namespace

// NOLINTNEXTLINE
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

  // NOLINTNEXTLINE
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
  for (const auto& d : vtk_file.point_data) {
    boost::apply_visitor(
        [&](auto e) {
          if (e.data.size() != vtk_file.unstructured_grid.points.size()) {
            throw base::LfException("Length of dataset " + e.name +
                                    " does not match the number of points.");
          }
        },
        d);
  }
  for (const auto& d : vtk_file.cell_data) {
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

/// Number of auxilliary nodes on this ref_el type (for order > 1)
unsigned int NumAuxNodes(base::RefEl ref_el, unsigned char order) {
  switch (ref_el) {
    case base::RefEl::kPoint():
      return 1;
    case base::RefEl::kSegment():
      return order - 1;
    case base::RefEl::kTria():
      return (2 - 3 * order + order * order) / 2;
    case base::RefEl::kQuad():
      return (order - 1) * (order - 1);
    default:
      LF_VERIFY_MSG(false, "RefElType " << ref_el << " not supported.");
  }
}

// Compute aux reference coordinates of lagrange points on segment
Eigen::MatrixXd AuxNodesSegment(unsigned char order) {
  LF_ASSERT_MSG(order > 1, "order <= 1");
  return Eigen::VectorXd::LinSpaced(order - 1, 1. / (order),
                                    (order - 1.) / order)
      .transpose();
}

// compute aux reference coordinates of lagrange points on Tria:
// NOLINTNEXTLINE(misc-no-recursion)
Eigen::MatrixXd AuxNodesTria(unsigned char order) {
  if (order < 3) {
    return {};
  }
  Eigen::MatrixXd result(2, NumAuxNodes(base::RefEl::kTria(), order));
  if (order == 3) {
    result << 1. / 3., 1. / 3.;
  } else if (order == 4) {
    result << 0.25, 0.5, 0.25, 0.25, 0.25, 0.5;
  } else {
    // assign nodes:
    result(0, 0) = 1. / order;
    result(1, 0) = 1. / order;
    result(0, 1) = (order - 2.) / order;
    result(1, 1) = 1. / order;
    result(0, 2) = 1. / order;
    result(1, 2) = (order - 2.) / order;

    if (order > 4) {
      auto segment_points = AuxNodesSegment(order - 3).eval();
      // assign edges:
      result.block(0, 3, 2, order - 4) =
          Eigen::Vector2d::UnitX() * segment_points * (order - 3.) / order +
          Eigen::Vector2d(1. / order, 1. / order) *
              Eigen::MatrixXd::Ones(1, order - 4);
      result.block(0, order - 1, 2, order - 4) =
          Eigen::Vector2d(-1, 1) * (order - 3.) / order * segment_points +
          Eigen::Vector2d((order - 2.) / order, 1. / order) *
              Eigen::MatrixXd::Ones(1, order - 4);
      result.block(0, 2 * order - 5, 2, order - 4) =
          -Eigen::Vector2d::UnitY() * (order - 3.) / order * segment_points +
          Eigen::Vector2d(1. / order, (order - 2.) / order) *
              Eigen::MatrixXd::Ones(1, order - 4);
    }
    if (order > 5) {
      // assign interior points recursively:
      auto points = AuxNodesTria(order - 3);
      result.block(0, 3 * order - 9, 2, points.cols()) =
          AuxNodesTria(order - 3) * (order - 6.) / order +
          Eigen::Vector2d(2. / order, 2. / order) *
              Eigen::MatrixXd::Ones(1, points.cols());
    }
  }
  return result;
}

// compute aux reference coordinates of lagrange points on quad:
Eigen::MatrixXd AuxNodesQuad(unsigned char order) {
  Eigen::MatrixXd result(2, NumAuxNodes(base::RefEl::kQuad(), order));
  result.row(0) =
      Eigen::VectorXd::LinSpaced(order - 1, 1. / order, (order - 1.) / order)
          .transpose()
          .replicate(1, order - 1);
  result.row(1) =
      Eigen::kroneckerProduct(Eigen::VectorXd::LinSpaced(order - 1, 1. / order,
                                                         (order - 1.) / order),
                              Eigen::VectorXd::Constant(order - 1, 1.))
          .transpose();
  return result;
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
VtkWriter::VtkWriter(std::shared_ptr<const mesh::Mesh> mesh,
                     std::string filename, dim_t codim, unsigned char order)
    : mesh_(std::move(mesh)),
      filename_(std::move(filename)),
      codim_(codim),
      order_(order),
      aux_node_offset_{{{mesh_, 0}, {mesh_, 1}}} {
  auto dim_mesh = mesh_->DimMesh();
  auto dim_world = mesh_->DimWorld();
  LF_ASSERT_MSG(dim_world > 0 && dim_world <= 4,
                "VtkWriter supports only dim_world = 1,2 or 3");
  LF_ASSERT_MSG(codim >= 0 && codim < dim_mesh, "codim out of bounds.");
  LF_ASSERT_MSG(order > 0 && order <= 10, "order must lie in [1,10]");

  // initialize reference coordinates of auxilliary nodes (if order_ > 0)
  if (order_ > 1) {
    aux_nodes_[base::RefEl::kSegment().Id()] = AuxNodesSegment(order);
    aux_nodes_[base::RefEl::kTria().Id()] = AuxNodesTria(order);
    aux_nodes_[base::RefEl::kQuad().Id()] = AuxNodesQuad(order);
  }

  // calculate total number of nodes (main + auxilliary nodes)
  unsigned int numNodes = mesh_->NumEntities(dim_mesh);
  for (auto ref_el :
       {base::RefEl::kSegment(), base::RefEl::kTria(), base::RefEl::kQuad()}) {
    if (ref_el.Dimension() > dim_mesh - codim_) {
      continue;
    }
    numNodes += NumAuxNodes(ref_el, order_) * mesh_->NumEntities(ref_el);
  }

  // insert main nodes:
  vtk_file_.unstructured_grid.points.resize(numNodes);
  Eigen::Matrix<double, 0, 1> zero;
  for (const auto* p : mesh_->Entities(dim_mesh)) {
    auto index = mesh_->Index(*p);
    Eigen::Vector3f coord;
    if (dim_world == 1) {
      coord(0) = static_cast<float>(p->Geometry()->Global(zero)(0));
      coord(1) = 0.F;
      coord(2) = 0.F;
    } else if (dim_world == 2) {
      coord.topRows<2>() = p->Geometry()->Global(zero).cast<float>();
      coord(2) = 0.F;
    } else {
      coord = p->Geometry()->Global(zero).cast<float>();
    }

    vtk_file_.unstructured_grid.points[index] = std::move(coord);
  }

  // insert auxilliary nodes:
  if (order > 1) {
    auto index_offset = mesh_->NumEntities(dim_mesh);
    for (char cd = static_cast<char>(dim_mesh - 1);
         cd >= static_cast<char>(codim); --cd) {
      for (const auto* e : mesh_->Entities(cd)) {
        auto ref_el = e->RefEl();
        if (ref_el == base::RefEl::kTria() && order < 3) {
          continue;
        }

        Eigen::MatrixXf coords(3, NumAuxNodes(ref_el, order));

        if (dim_world == 1) {
          coords.row(0) =
              e->Geometry()->Global(aux_nodes_[ref_el.Id()]).cast<float>();
          coords.row(1).setZero();
          coords.row(2).setZero();
        } else if (dim_world == 2) {
          coords.topRows(2) =
              e->Geometry()->Global(aux_nodes_[ref_el.Id()]).cast<float>();
          coords.row(2).setZero();
        } else {
          coords = e->Geometry()->Global(aux_nodes_[ref_el.Id()]).cast<float>();
        }
        for (Eigen::Index i = 0; i < coords.cols(); ++i) {
          vtk_file_.unstructured_grid.points[index_offset + i] = coords.col(i);
        }
        aux_node_offset_[cd](*e) = index_offset;
        index_offset += coords.cols();
      }
    }
  }

  // compute number of main/aux nodes for each element:
  std::array<unsigned int, 5> num_nodes{};
  num_nodes[base::RefEl::kSegment().Id()] =
      2 + NumAuxNodes(base::RefEl::kSegment(), order);
  num_nodes[base::RefEl::kTria().Id()] =
      3 + 3 * NumAuxNodes(base::RefEl::kSegment(), order) +
      NumAuxNodes(base::RefEl::kTria(), order);
  num_nodes[base::RefEl::kQuad().Id()] =
      4 + 4 * NumAuxNodes(base::RefEl::kSegment(), order) +
      NumAuxNodes(base::RefEl::kQuad(), order);

  // insert elements:
  vtk_file_.unstructured_grid.cells.resize(mesh_->NumEntities(codim));
  vtk_file_.unstructured_grid.cell_types.resize(mesh_->NumEntities(codim));
  auto points_per_segment = NumAuxNodes(base::RefEl::kSegment(), order);
  for (const auto* e : mesh_->Entities(codim)) {
    auto index = mesh_->Index(*e);
    auto ref_el = e->RefEl();
    auto& node_indices = vtk_file_.unstructured_grid.cells[index];
    node_indices.reserve(num_nodes[ref_el.Id()]);

    // node indices that make up this cell:
    for (const auto* p : e->SubEntities(dim_mesh - codim)) {
      node_indices.push_back(mesh_->Index(*p));
    }

    // node indices of segments of this cell:
    auto addSegmentNodes = [&](const mesh::Entity& s, bool invert) {
      auto start_index = aux_node_offset_[dim_mesh - 1](s);
      if (!invert) {
        for (unsigned int i = 0; i < points_per_segment; ++i) {
          node_indices.push_back(start_index + i);
        }
      } else {
        for (int i = static_cast<int>(points_per_segment - 1); i >= 0; --i) {
          node_indices.push_back(start_index + i);
        }
      }
    };
    switch (ref_el) {
      case base::RefEl::kSegment():
        addSegmentNodes(*e, false);
        break;
      case base::RefEl::kTria(): {
        const auto* iterator = e->SubEntities(1).begin();
        const auto* o_iterator = e->RelativeOrientations().begin();
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::negative);
        ++iterator;
        ++o_iterator;
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::negative);
        ++iterator;
        ++o_iterator;
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::negative);
        break;
      }
      case base::RefEl::kQuad(): {
        const auto* iterator = e->SubEntities(1).begin();
        const auto* o_iterator = e->RelativeOrientations().begin();
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::negative);
        ++iterator;
        ++o_iterator;
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::negative);
        ++iterator;
        ++o_iterator;
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::positive);
        ++iterator;
        ++o_iterator;
        addSegmentNodes(**iterator,
                        (*o_iterator) == mesh::Orientation::positive);
        break;
      }
      default:
        LF_VERIFY_MSG(false, "something is wrong");
    }

    // indices in the interior of the cell:
    if (dim_mesh - codim > 1) {
      auto offset = aux_node_offset_[0](*e);
      for (unsigned i = 0; i < NumAuxNodes(e->RefEl(), order); ++i) {
        node_indices.push_back(offset + i);
      }
    }

    if (order_ == 1) {
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
    } else {
      switch (ref_el) {
        case base::RefEl::kSegment():
          vtk_file_.unstructured_grid.cell_types[index] =
              VtkFile::CellType::VTK_LAGRANGE_CURVE;
          break;
        case base::RefEl::kTria():
          vtk_file_.unstructured_grid.cell_types[index] =
              VtkFile::CellType::VTK_LAGRANGE_TRIANGLE;
          break;
        case base::RefEl::kQuad():
          vtk_file_.unstructured_grid.cell_types[index] =
              VtkFile::CellType::VTK_LAGRANGE_QUADRILATERAL;
          break;
        default:
          LF_VERIFY_MSG(false, "Unsupported RefElType " << ref_el);
      }
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

void VtkWriter::WritePointData(
    const std::string& name, const mesh::utils::MeshDataSet<unsigned int>& mds,
    unsigned int undefined_value) {
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

void VtkWriter::WritePointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::VectorXd>& mds,
    const Eigen::VectorXd& undefined_value) {
  WriteVectorPointData<Eigen::Dynamic, double>(name, mds, undefined_value);
}

void VtkWriter::WritePointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::VectorXf>& mds,
    const Eigen::VectorXf& undefined_value) {
  WriteVectorPointData<Eigen::Dynamic, float>(name, mds, undefined_value);
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
                              const mesh::utils::MeshDataSet<unsigned int>& mds,
                              unsigned int undefined_value) {
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

void VtkWriter::WriteCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::VectorXd>& mds,
    const Eigen::VectorXd& undefined_value) {
  WriteVectorCellData<Eigen::Dynamic, double>(name, mds, undefined_value);
}

void VtkWriter::WriteCellData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::VectorXf>& mds,
    const Eigen::VectorXf& undefined_value) {
  WriteVectorCellData<Eigen::Dynamic, float>(name, mds, undefined_value);
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

template <class T>
void VtkWriter::WriteScalarPointData(const std::string& name,
                                     const mesh::utils::MeshDataSet<T>& mds,
                                     T undefined_value) {
  LF_ASSERT_MSG(order_ == 1,
                "WritePointData accepts MeshDataSets only if order = 1. For "
                "order > 1 you have to provide MeshFunctions.");
  CheckAttributeSetName(vtk_file_.point_data, name);
  VtkFile::ScalarData<T> data{};
  data.data.resize(mesh_->NumEntities(mesh_->DimMesh()));
  data.name = name;
  for (const auto* p : mesh_->Entities(mesh_->DimMesh())) {
    if (mds.DefinedOn(*p)) {
      data.data[mesh_->Index(*p)] = mds(*p);
    } else {
      data.data[mesh_->Index(*p)] = undefined_value;
    }
  }
  vtk_file_.point_data.push_back(std::move(data));
}

template <int ROWS, class T>
void VtkWriter::WriteVectorPointData(
    const std::string& name,
    const mesh::utils::MeshDataSet<Eigen::Matrix<T, ROWS, 1>>& mds,
    const Eigen::Matrix<T, ROWS, 1>& undefined_value) {
  LF_ASSERT_MSG(order_ == 1,
                "WritePointData accepts MeshDataSets only if order = 1. For "
                "order > 1 you have to provide MeshFunctions.");
  CheckAttributeSetName(vtk_file_.point_data, name);
  VtkFile::VectorData<T> data{};
  data.data.resize(mesh_->NumEntities(mesh_->DimMesh()));
  data.name = name;
  Eigen::Matrix<T, 3, 1> undefined_value_padded;
  PadWithZeros<ROWS, T>(undefined_value_padded, undefined_value);

  for (const auto* p : mesh_->Entities(mesh_->DimMesh())) {
    if (mds.DefinedOn(*p)) {
      PadWithZeros<ROWS, T>(data.data[mesh_->Index(*p)], mds(*p));
    } else {
      data.data[mesh_->Index(*p)] = undefined_value_padded;
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
  data.data.resize(mesh_->NumEntities(codim_));
  data.name = name;
  for (const auto* e : mesh_->Entities(codim_)) {
    if (mds.DefinedOn(*e)) {
      data.data[mesh_->Index(*e)] = mds(*e);
    } else {
      data.data[mesh_->Index(*e)] = undefined_value;
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
  data.data.resize(mesh_->NumEntities(codim_));
  data.name = name;
  Eigen::Matrix<T, 3, 1> undefined_value_padded;
  PadWithZeros<ROWS, T>(undefined_value_padded, undefined_value);

  for (const auto* p : mesh_->Entities(codim_)) {
    if (mds.DefinedOn(*p)) {
      PadWithZeros<ROWS, T>(data.data[mesh_->Index(*p)], mds(*p));
    } else {
      data.data[mesh_->Index(*p)] = undefined_value_padded;
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
