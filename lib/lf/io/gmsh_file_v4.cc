/**
 * @file
 * @brief Implementation of all logic related to the class GmshFileV4
 * @author Raffael Casagrande
 * @date   2019-09-02 03:02:19
 * @copyright MIT License
 */

#include "gmsh_file_v4.h"

namespace lf::io {

/// Output the element type onto a stream:
std::ostream& operator<<(std::ostream& stream, GMshFileV4::ElementType et) {
  switch (et) {
    default:
      break;
    case GMshFileV4::ElementType::EDGE2:
      stream << "EDGE2";
      break;
    case GMshFileV4::ElementType::TRIA3:
      stream << "TRIA3";
      break;
    case GMshFileV4::ElementType::QUAD4:
      stream << "QUAD4";
      break;
    case GMshFileV4::ElementType::TET4:
      stream << "TET4";
      break;
    case GMshFileV4::ElementType::HEX8:
      stream << "HEX8";
      break;
    case GMshFileV4::ElementType::PRISM6:
      stream << "PRISM6";
      break;
    case GMshFileV4::ElementType::PYRAMID5:
      stream << "PYRAMID5";
      break;
    case GMshFileV4::ElementType::EDGE3:
      stream << "EDGE3";
      break;
    case GMshFileV4::ElementType::TRIA6:
      stream << "TRIA6";
      break;
    case GMshFileV4::ElementType::QUAD9:
      stream << "QUAD9";
      break;
    case GMshFileV4::ElementType::TET10:
      stream << "TET10";
      break;
    case GMshFileV4::ElementType::HEX27:
      stream << "HEX27";
      break;
    case GMshFileV4::ElementType::PRISM18:
      stream << "PRISM18";
      break;
    case GMshFileV4::ElementType::PYRAMID14:
      stream << "PYRAMID14";
      break;
    case GMshFileV4::ElementType::POINT:
      stream << "POINT";
      break;
    case GMshFileV4::ElementType::QUAD8:
      stream << "QUAD8";
      break;
    case GMshFileV4::ElementType::HEX20:
      stream << "HEX20";
      break;
    case GMshFileV4::ElementType::PRISM15:
      stream << "PRISM15";
      break;
    case GMshFileV4::ElementType::PYRAMID13:
      stream << "PYRAMID13";
      break;
    case GMshFileV4::ElementType::TRIA9:
      stream << "TRIA9";
      break;
    case GMshFileV4::ElementType::TRIA10:
      stream << "TRIA10";
      break;
    case GMshFileV4::ElementType::TRIA12:
      stream << "TRIA12";
      break;
    case GMshFileV4::ElementType::TRIA15:
      stream << "TRIA15";
      break;
    case GMshFileV4::ElementType::TRIA15_5:
      stream << "TRIA15_5";
      break;
    case GMshFileV4::ElementType::TRIA21:
      stream << "TRIA21";
      break;
    case GMshFileV4::ElementType::EDGE4:
      stream << "EDGE4";
      break;
    case GMshFileV4::ElementType::EDGE5:
      stream << "EDGE5";
      break;
    case GMshFileV4::ElementType::EDGE6:
      stream << "EDGE6";
      break;
    case GMshFileV4::ElementType::TET20:
      stream << "TET20";
      break;
    case GMshFileV4::ElementType::TET35:
      stream << "TET35";
      break;
    case GMshFileV4::ElementType::TET56:
      stream << "TET56";
      break;
    case GMshFileV4::ElementType::HEX64:
      stream << "HEX64";
      break;
    case GMshFileV4::ElementType::HEX125:
      stream << "HEX125";
      break;
  }
  return stream;
}

/// Number of nodes that this element type has
int NumNodes(GMshFileV4::ElementType et) {
  switch (et) {
    default:
      break;
    case GMshFileV4::ElementType::EDGE2:
      return 2;
    case GMshFileV4::ElementType::TRIA3:
      return 3;
    case GMshFileV4::ElementType::QUAD4:
    case GMshFileV4::ElementType::TET4:
      return 4;
    case GMshFileV4::ElementType::HEX8:
      return 8;
    case GMshFileV4::ElementType::PRISM6:
      return 6;
    case GMshFileV4::ElementType::PYRAMID5:
      return 5;
    case GMshFileV4::ElementType::EDGE3:
      return 3;
    case GMshFileV4::ElementType::TRIA6:
      return 6;
    case GMshFileV4::ElementType::QUAD9:
      return 9;
    case GMshFileV4::ElementType::TET10:
      return 10;
    case GMshFileV4::ElementType::HEX27:
      return 27;
    case GMshFileV4::ElementType::PRISM18:
      return 18;
    case GMshFileV4::ElementType::PYRAMID14:
      return 14;
    case GMshFileV4::ElementType::POINT:
      return 1;
    case GMshFileV4::ElementType::QUAD8:
      return 8;
    case GMshFileV4::ElementType::HEX20:
      return 20;
    case GMshFileV4::ElementType::PRISM15:
      return 15;
    case GMshFileV4::ElementType::PYRAMID13:
      return 13;
    case GMshFileV4::ElementType::TRIA9:
      return 9;
    case GMshFileV4::ElementType::TRIA10:
      return 10;
    case GMshFileV4::ElementType::TRIA12:
      return 12;
    case GMshFileV4::ElementType::TRIA15:
    case GMshFileV4::ElementType::TRIA15_5:
      return 15;
    case GMshFileV4::ElementType::TRIA21:
      return 21;
    case GMshFileV4::ElementType::EDGE4:
      return 4;
    case GMshFileV4::ElementType::EDGE5:
      return 5;
    case GMshFileV4::ElementType::EDGE6:
      return 6;
    case GMshFileV4::ElementType::TET20:
      return 20;
    case GMshFileV4::ElementType::TET35:
      return 35;
    case GMshFileV4::ElementType::TET56:
      return 56;
    case GMshFileV4::ElementType::HEX64:
      return 64;
    case GMshFileV4::ElementType::HEX125:
      return 125;
  }
  LF_VERIFY_MSG(false, "unknown Gmsh element type");
  // Make compiler happy:
  return 0;
}

base::RefEl RefElOf(GMshFileV4::ElementType et) {
  switch (et) {
    case GMshFileV4::ElementType::POINT:
      return base::RefEl::kPoint();

    case GMshFileV4::ElementType::EDGE2:
    case GMshFileV4::ElementType::EDGE3:
    case GMshFileV4::ElementType::EDGE4:
    case GMshFileV4::ElementType::EDGE5:
    case GMshFileV4::ElementType::EDGE6:
      return base::RefEl::kSegment();

    case GMshFileV4::ElementType::TRIA3:
    case GMshFileV4::ElementType::TRIA6:
    case GMshFileV4::ElementType::TRIA9:
    case GMshFileV4::ElementType::TRIA10:
    case GMshFileV4::ElementType::TRIA12:
    case GMshFileV4::ElementType::TRIA15:
    case GMshFileV4::ElementType::TRIA15_5:
    case GMshFileV4::ElementType::TRIA21:
      return base::RefEl::kTria();

    case GMshFileV4::ElementType::QUAD4:
    case GMshFileV4::ElementType::QUAD8:
    case GMshFileV4::ElementType::QUAD9:
      return base::RefEl::kQuad();

    case GMshFileV4::ElementType::TET4:
    case GMshFileV4::ElementType::HEX8:
    case GMshFileV4::ElementType::PRISM6:
    case GMshFileV4::ElementType::PYRAMID5:
    case GMshFileV4::ElementType::TET10:
    case GMshFileV4::ElementType::HEX27:
    case GMshFileV4::ElementType::PRISM18:
    case GMshFileV4::ElementType::PYRAMID14:
    case GMshFileV4::ElementType::HEX20:
    case GMshFileV4::ElementType::PRISM15:
    case GMshFileV4::ElementType::PYRAMID13:
    case GMshFileV4::ElementType::TET20:
    case GMshFileV4::ElementType::TET35:
    case GMshFileV4::ElementType::TET56:
    case GMshFileV4::ElementType::HEX64:
    case GMshFileV4::ElementType::HEX125:
    default:
      LF_VERIFY_MSG(
          false, "Reference element not supported for GmshElement type " << et);
  }
}

/// Dimension of the GmshElement type
int DimOf(GMshFileV4::ElementType et) {
  switch (et) {
    case GMshFileV4::ElementType::POINT:
      return 0;
    case GMshFileV4::ElementType::EDGE2:
    case GMshFileV4::ElementType::EDGE3:
    case GMshFileV4::ElementType::EDGE4:
    case GMshFileV4::ElementType::EDGE5:
    case GMshFileV4::ElementType::EDGE6:
      return 1;
    case GMshFileV4::ElementType::TRIA3:
    case GMshFileV4::ElementType::QUAD4:
    case GMshFileV4::ElementType::TRIA6:
    case GMshFileV4::ElementType::QUAD9:
    case GMshFileV4::ElementType::QUAD8:
    case GMshFileV4::ElementType::TRIA9:
    case GMshFileV4::ElementType::TRIA10:
    case GMshFileV4::ElementType::TRIA12:
    case GMshFileV4::ElementType::TRIA15:
    case GMshFileV4::ElementType::TRIA15_5:
    case GMshFileV4::ElementType::TRIA21:
      return 2;
    case GMshFileV4::ElementType::TET4:
    case GMshFileV4::ElementType::HEX8:
    case GMshFileV4::ElementType::PRISM6:
    case GMshFileV4::ElementType::PYRAMID5:
    case GMshFileV4::ElementType::TET10:
    case GMshFileV4::ElementType::HEX27:
    case GMshFileV4::ElementType::PRISM18:
    case GMshFileV4::ElementType::PYRAMID14:
    case GMshFileV4::ElementType::HEX20:
    case GMshFileV4::ElementType::PRISM15:
    case GMshFileV4::ElementType::PYRAMID13:
    case GMshFileV4::ElementType::TET20:
    case GMshFileV4::ElementType::TET35:
    case GMshFileV4::ElementType::TET56:
    case GMshFileV4::ElementType::HEX64:
    case GMshFileV4::ElementType::HEX125:
      return 3;
    default:
      LF_VERIFY_MSG(false, "Unknown GmshElement Type.");
  }
  // Make compiler happy:
  return -1;
}

namespace detail {

// defined in gmsh_file_v4_text.cc
bool ParseGmshFileV4Text(std::string::const_iterator begin,
                         std::string::const_iterator end, GMshFileV4* result);

bool ParseGmshFileV4Binary(std::string::const_iterator begin,
                           std::string::const_iterator end, int one,
                           GMshFileV4* result);

}  // namespace detail

GMshFileV4 ReadGmshFileV4(std::string::const_iterator begin,
                          std::string::const_iterator end,
                          const std::string& version, bool is_binary,
                          int size_t_size, int one,
                          const std::string& filename) {
  LF_VERIFY_MSG(version == "4.1", "Only version 4.1 is supported so far");
  LF_VERIFY_MSG(size_t_size == sizeof(std::size_t),
                "size of size_t must be " << sizeof(std::size_t));

  GMshFileV4 result;
  result.version_number = version;
  result.is_binary = is_binary;
  result.size_t_size = size_t_size;

  bool succesful = false;
  if (!is_binary) {
    succesful = detail::ParseGmshFileV4Text(begin, end, &result);
  } else {
    succesful = detail::ParseGmshFileV4Binary(begin, end, one, &result);
  }

  LF_VERIFY_MSG(succesful, "Could not parse file " << filename);
  // LF_VERIFY_MSG(iter == end, "Could not parse all of file " << filename);

  // transpose all periodic matrices because they are read in column-first mode
  // but gmsh writes them in row-first mode
  for (auto& p : result.periodic_links) {
    if (p.affine_transform) {
      p.affine_transform->transposeInPlace();
    }
  }

  return result;
}

}  // namespace lf::io
