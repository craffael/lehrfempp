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

  VariantVisitor(LAMBDA lambda) : lambda_(lambda) {}

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
void cell_list_size(unsigned int& result,
                    const std::vector<std::vector<unsigned int>>& cells) {
  result = cells.size();
  for (auto& v : cells) {
    result += v.size();
  }
}
}  // namespace

BOOST_PHOENIX_ADAPT_FUNCTION(void, CellListSize, cell_list_size, 2);

// Implement the container_size custom generator
////////////////////////////////////////////////////////////////////////////////

// namespace lf::io {
// BOOST_SPIRIT_TERMINAL(container_size);
// }
//
// namespace boost::spirit {
// template <>
// struct use_terminal<boost::spirit::karma::domain,
// lf::io::tag::container_size>
//     : boost::mpl::true_ {};
// }  // namespace boost::spirit
//
// namespace lf::io {
// template <typename Subject>
// struct container_size_generator : boost::spirit::karma::primitive_generator<
//                                       container_size_generator<Subject>> {
//   using properties = typename Subject::properties;
//   Subject subject;
//
//   container_size_generator(const Subject& s) : subject(s) {}
//
//   template <class Context, class Iterator>
//   struct attribute {
//     using type = typename Subject::size_type;
//   };
//
//   template<class Iterator, class Context, class Delimiter, class
//   Attribute> bool generate(Iterator& sink, Context& ctx, const
//   Delimiter& delimiter, const Attribute& attr) const {
//     return boost::spirit::karma::int_inserter<10, >
//   }
// };
// }  // namespace lf::io

// End of container_size custom generator
////////////////////////////////////////////////////////////////////////////////

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

    if constexpr (BINARY) {
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
    if constexpr (BINARY) {
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

    if constexpr (BINARY) {
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
    if constexpr (BINARY) {
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
    if constexpr (BINARY) {
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
}  // namespace lf::io
