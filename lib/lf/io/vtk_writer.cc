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

// clang-format off
BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile::UnstructuredGrid,
    (std::vector<Eigen::Vector3d>, points)
    (std::vector<std::vector<unsigned int>>, cells)
    (std::vector<lf::io::VtkFile::CellType>, cell_types)
);

BOOST_FUSION_ADAPT_TPL_STRUCT((T), (lf::io::VtkFile::ScalarData)(T),
    (std::string, data_name)
    (std::string, lookup_table)
    (auto, data)
);

 BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile::Attributes,
     (std::vector<lf::io::VtkFile::ScalarData<int>>, scalar_int_data)
     (std::vector<lf::io::VtkFile::ScalarData<unsigned int>>, scalar_unsigned_int_data)
     (std::vector<lf::io::VtkFile::ScalarData<long>>, scalar_long_data)
     (std::vector<lf::io::VtkFile::ScalarData<unsigned long>>, scalar_unsigned_long_data)
     (std::vector<lf::io::VtkFile::ScalarData<float>>, scalar_float_data)
     (std::vector<lf::io::VtkFile::ScalarData<double>>, scalar_double_data)
 );

BOOST_FUSION_ADAPT_STRUCT(lf::io::VtkFile,
    (std::string, header)
    (lf::io::VtkFile::Format, format)
    (lf::io::VtkFile::UnstructuredGrid, unstructured_grid)
    (lf::io::VtkFile::Attributes, point_data)
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
//   template<class OutputIterator, class Context, class Delimiter, class
//   Attribute> bool generate(OutputIterator& sink, Context& ctx, const
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
template <class OutputIterator, bool BINARY>
struct VtkGrammar : karma::grammar<OutputIterator, VtkFile()> {
  // template <class T>
  // karma::rule<OutputIterator, std::vector<T>()> vector_size() {
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
      points_vector %= *(karma::big_bin_double << karma::big_bin_double
                                               << karma::big_bin_double)
                       << eol;
    } else {
      points_vector %= *(karma::double_ << ' ' << karma::double_ << ' '
                                        << karma::double_ << karma::eol);
    }

    points = lit("POINTS ")
             << karma::uint_[karma::_1 = boost::phoenix::size(karma::_val)]
             << lit(" double") << karma::eol
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

    scalar_data_int %= (lit("SCALARS ") << karma::string << lit(" int 1") << eol
                                        << lit("LOOKUP_TABLE ") << karma::string
                                        << eol << (*(karma::int_ % eol)));
    scalar_data_unsigned_int %=
        (lit("SCALARS ") << karma::string << lit(" unsigned_int 1") << eol
                         << lit("LOOKUP_TABLE ") << karma::string << eol
                         << (*(karma::uint_ % eol)));

    scalar_data_long %=
        (lit("SCALARS ") << karma::string << lit(" long 1") << eol
                         << lit("LOOKUP_TABLE ") << karma::string << eol
                         << (*(karma::long_ % eol)));

    scalar_data_unsigned_long %=
        (lit("SCALARS ") << karma::string << lit(" unsigned_long 1") << eol
                         << lit("LOOKUP_TABLE ") << karma::string << eol
                         << (*(karma::ulong_ % eol)));

    scalar_data_float %=
        (lit("SCALARS ") << karma::string << lit(" float 1") << eol
                         << lit("LOOKUP_TABLE ") << karma::string << eol
                         << (*(karma::float_ % eol)));

    scalar_data_double %=
        (lit("SCALARS ") << karma::string << lit(" double 1") << eol
                         << lit("LOOKUP_TABLE ") << karma::string << eol
                         << (*(karma::double_ % eol)));

    attributes %= (*scalar_data_int)
                  << (*scalar_data_unsigned_int) << (*scalar_data_long)
                  << (*scalar_data_unsigned_long) << (*scalar_data_float)
                  << (*scalar_data_double);

    root %= lit("# vtk DataFile Version 3.0")
            << karma::eol << karma::string << karma::eol << karma::stream
            << karma::eol << unstructured_grid << karma::eol;
    // << lit("POINT_DATA ") << lit(size(at_c<0>(at_c<2>(_val)))) << eol
    // << attributes;
  }

  karma::rule<OutputIterator, VtkFile()> root;
  karma::rule<OutputIterator, VtkFile()> header;
  karma::rule<OutputIterator, VtkFile::UnstructuredGrid()> unstructured_grid;
  karma::rule<OutputIterator, std::vector<Eigen::Vector3d>()> points;
  karma::rule<OutputIterator, std::vector<Eigen::Vector3d>()> points_vector;
  karma::rule<OutputIterator, std::vector<std::vector<unsigned int>>> cells;
  karma::rule<OutputIterator, std::vector<unsigned int>> cells_vector;
  karma::rule<OutputIterator, std::vector<VtkFile::CellType>> cell_types;
  karma::rule<OutputIterator, VtkFile::CellType> cell_type;
  karma::rule<OutputIterator, VtkFile::Attributes()> attributes;
  karma::rule<OutputIterator, VtkFile::ScalarData<int>()> scalar_data_int;
  karma::rule<OutputIterator, VtkFile::ScalarData<unsigned int>()>
      scalar_data_unsigned_int;
  karma::rule<OutputIterator, VtkFile::ScalarData<long>()> scalar_data_long;
  karma::rule<OutputIterator, VtkFile::ScalarData<unsigned long>()>
      scalar_data_unsigned_long;
  karma::rule<OutputIterator, VtkFile::ScalarData<float>()> scalar_data_float;
  karma::rule<OutputIterator, VtkFile::ScalarData<double>()> scalar_data_double;
  karma::rule<OutputIterator, std::vector<int>()> temp;
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

  VtkGrammar<decltype(outit), false> grammar{};
  bool result = karma::generate(outit, grammar, vtk_file);
  file.close();
  if (!result) {
    throw base::LfException("Karma error");
  }
}
}  // namespace lf::io
