#include "gmsh_reader.h"

#include <lf/geometry/geometry.h>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_string.hpp>
#include "eigen_fusion_adapter.h"

#include <fstream>

using size_type = lf::mesh::Mesh::size_type;

// Structures that represent the MshFile:
namespace lf::io {

bool GmshReader::IsPhysicalEntity(const mesh::Entity& e,
                                  size_type physical_entity_nr) const {
  auto physical_entities = PhysicalEntityNr(e);
  return std::find(physical_entities.begin(), physical_entities.end(),
                   physical_entity_nr) != physical_entities.end();
}

GmshReader::GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
                       const ::lf::io::GMshFileV2& msh_file)
    : mesh_factory_(std::move(factory)) {
  // 1) Check Gmsh_file and initialize
  //////////////////////////////////////////////////////////////////////////

  dim_t dim_mesh = mesh_factory_->DimMesh();
  dim_t dim_world = mesh_factory_->DimWorld();
  LF_VERIFY_MSG(
      dim_mesh >= 2 && dim_mesh <= 3 && dim_world >= 2 && dim_world <= 3,
      "GmshReader supports only 2D and 3D meshes.");

  // 1) Determine which nodes of gmsh are also nodes of the LehrFEM++ mesh +
  // count number of entities of each codimension (exclude auxilliary nodes).
  // This is necessary because e.g. second order meshes in gmsh need a lot of
  // auxilliary nodes to define their geometry...
  /////////////////////////////////////////////////////////////////////////////
  // is_main_node[i] = true means that gmsh-node i is one of the main nodes of
  // an element
  std::vector<bool> is_main_node(msh_file.Nodes.size(), false);

  // mi2gi = mesh_index_2_gmsh_index
  // mi2gi[c][i] contains the gmsh entities that belong to the mesh entity with
  //             codim = c and mesh index i.
  std::vector<std::vector<std::vector<size_type>>> mi2gi(dim_mesh + 1);

  {
    // count the number of entities for each codimension and reserve space:
    std::vector<size_type> num_entities(mesh_factory_->DimMesh() + 1, 0);

    for (const auto& e : msh_file.Elements) {
      LF_ASSERT_MSG(DimOf(e.Type) <= dim_mesh,
                    "mesh_factory->DimMesh() = "
                        << dim_mesh
                        << ", but msh-file contains entities with dimension "
                        << DimOf(e.Type));

      ++num_entities[DimOf(e.Type)];

      if (DimOf(e.Type) == dim_mesh) {
        // mark main nodes
        auto ref_el = RefElOf(e.Type);
        for (unsigned int i = 0; i < ref_el.NumNodes(); ++i) {
          auto node_number = e.NodeNumbers[i];
          if (is_main_node.size() <= node_number) {
            is_main_node.resize(node_number + 1);
          }
          is_main_node[e.NodeNumbers[i]] = true;
        }
      }
    }

    for (dim_t c = 0; c <= dim_mesh; ++c) {
      mi2gi[c].reserve(num_entities[dim_mesh - c]);
    }

    LF_ASSERT_MSG(num_entities[dim_mesh] > 0,
                  "MshFile contains no elements with dimension " << dim_mesh);
  }

  // 2) Insert main nodes into MeshFactory
  /////////////////////////////////////////////////////////////////////////////

  // gmsh_index_2_mesh_index for nodes:
  // gi2mi[i] = j means that gmsh node with gmsh index i has mesh index j
  std::vector<size_type> gi2mi;
  gi2mi.resize(is_main_node.size(), -1);

  // gi2i[i] = j implies msh_file.Nodes[j].first = i
  std::vector<size_type> gi2i;
  gi2i.resize(is_main_node.size(), -1);

  for (std::size_t i = 0; i < msh_file.Nodes.size(); ++i) {
    auto& n = msh_file.Nodes[i];
    if (gi2i.size() <= n.first) {
      gi2i.resize(n.first + 1, -1);
    }
    gi2i[n.first] = i;

    if (is_main_node.size() <= n.first || !is_main_node[n.first]) {
      continue;
    }

    size_type mi;
    if (dim_world == 2) {
      LF_ASSERT_MSG(
          n.second(2) == 0,
          "In a 2D GmshMesh, the z-coordinate of every node must be zero");
      mi = mesh_factory_->AddPoint(n.second.topRows(2));
    } else {
      mi = mesh_factory_->AddPoint(n.second);
    }
    gi2mi[n.first] = mi;
  }

  // 3) Insert entities (except nodes) into MeshFactory:
  //////////////////////////////////////////////////////////////////////////////

  std::size_t begin = 0;
  for (std::size_t end = 0; end < msh_file.Elements.size(); ++end) {
    auto& begin_element = msh_file.Elements[begin];
    auto& end_element = msh_file.Elements[end];
    auto ref_el = RefElOf(end_element.Type);
    auto codim = dim_mesh - ref_el.Dimension();
    if (begin_element.NodeNumbers == end_element.NodeNumbers && begin != end &&
        begin_element.Type == end_element.Type) {
      // This entity appears more than once
      mi2gi[codim].back().push_back(end);
      continue;
    }

    begin = end;
    if (ref_el == base::RefEl::kPoint()) {
      // special case, this entity is a point (which has already been inserted)
      auto mesh_index = gi2mi[end_element.NodeNumbers[0]];
      if (mi2gi[dim_mesh].size() <= mesh_index) {
        mi2gi[dim_mesh].resize(mesh_index + 1);
      }
      mi2gi[dim_mesh][mesh_index].push_back(end);
    } else {
      // gmsh element is not a point -> insert entity:
      auto num_nodes = end_element.NodeNumbers.size();
      Eigen::MatrixXd node_coords(dim_world, num_nodes);
      for (std::size_t i = 0; i < num_nodes; ++i) {
        auto node_coord =
            msh_file.Nodes[gi2i[end_element.NodeNumbers[i]]].second;
        if (dim_world == 2) {
          node_coords.col(i) = node_coord.topRows(2);
        } else {
          node_coords.col(i) = node_coord;
        }
      }
      std::unique_ptr<geometry::Geometry> geom;

      switch (end_element.Type) {
        case GMshFileV2::ElementType::EDGE2:
          ref_el = base::RefEl::kSegment();
          geom = std::make_unique<geometry::SegmentO1>(node_coords);
          break;
        case GMshFileV2::ElementType::EDGE3:
          ref_el = base::RefEl::kSegment();
          geom = std::make_unique<geometry::SegmentO2>(node_coords);
          break;
        case GMshFileV2::ElementType::TRIA3:
          ref_el = base::RefEl::kTria();
          geom = std::make_unique<geometry::TriaO1>(node_coords);
          break;
        case GMshFileV2::ElementType::TRIA6:
          ref_el = base::RefEl::kTria();
          geom = std::make_unique<geometry::TriaO2>(node_coords);
          break;
        case GMshFileV2::ElementType::QUAD4:
          ref_el = base::RefEl::kQuad();
          geom = std::make_unique<geometry::QuadO1>(node_coords);
          break;
        case GMshFileV2::ElementType::QUAD8:
          ref_el = base::RefEl::kQuad();
          geom = std::make_unique<geometry::QuadO2>(node_coords);
          break;
        case GMshFileV2::ElementType::QUAD9:
          ref_el = base::RefEl::kQuad();
          geom = std::make_unique<geometry::QuadO2>(node_coords.leftCols(8));
          break;
        default:
          LF_VERIFY_MSG(false, "Gmsh element type "
                                   << end_element.Type
                                   << " not (yet) supported by GmshReader.");
      }
      std::vector<size_type> main_nodes(ref_el.NumNodes());
      for (dim_t i = 0; i < ref_el.NumNodes(); ++i) {
        main_nodes[i] = gi2mi[end_element.NodeNumbers[i]];
      }

      mesh_factory_->AddEntity(ref_el, main_nodes, std::move(geom));
      mi2gi[codim].emplace_back(std::vector{static_cast<unsigned int>(end)});
    }
  }

  // 4) Construct mesh
  //////////////////////////////////////////////////////////////////////////////
  mesh_ = mesh_factory_->Build();

  // 5) Build MeshDataSet that assigns the physical entitiies:
  //////////////////////////////////////////////////////////////////////////////
  physical_nrs_ =
      mesh::utils::make_AllCodimMeshDataSet<std::vector<size_type>>(mesh_);

  for (dim_t c = 0; c <= dim_mesh; ++c) {
    for (auto& e : mesh_->Entities(c)) {
      auto mi = mesh_->Index(e);
      if (c == dim_mesh && mi >= mi2gi[dim_mesh].size()) {
        // this point did not appear as a gmsh element in the file -> don't
        // assign any physical entity nr.
        continue;
      }
      if (mi2gi[c].size() > mi) {
        std::vector<size_type> temp;
        for (auto& gmsh_index : mi2gi[c][mi]) {
          temp.push_back(msh_file.Elements[gmsh_index].PhysicalEntityNr);
        }

        physical_nrs_->operator()(e) = std::move(temp);
      }
    }
  }

  // 6) Create mapping physicalEntityNr <-> physicalEntityName:
  //////////////////////////////////////////////////////////////////////////

  for (auto& pe : msh_file.PhysicalEntities) {
    name_2_nr_.insert(
        std::pair{pe.Name, std::pair{pe.Number, dim_mesh - pe.Dimension}});
    nr_2_name_.insert(
        std::pair{pe.Number, std::pair{pe.Name, dim_mesh - pe.Dimension}});
  }

  if (!msh_file.Periodic.empty()) {
    /*LOGGER_ENTRY(logger_,
                 "WARNING: GMSH File  contains periodic boundary relations "
                 "between elements. These are ignored by GmshReader.",
                 3);*/
  }
}

GmshReader::GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
                       const std::string& filename)
    : GmshReader(std::move(factory), readGmshFileV2(filename)) {}

size_type GmshReader::PhysicalEntityName2Nr(const std::string& name,
                                            dim_t codim) const {
  LF_ASSERT_MSG(!name.empty(), "name is empty");
  auto [begin, end] = name_2_nr_.equal_range(name);  // NOLINT
  if (begin == end) {
    throw base::LfException("No Physical Entity with this name found.");
  }
  auto result = *begin;
  ++begin;
  if (begin == end) {
    if (codim == dim_t(-1) || codim == result.second.second) {
      return result.second.first;
    }
  } else {
    if (codim == dim_t(-1)) {
      throw base::LfException(
          "There are multiple physical entities with the name " + name +
          ", please specify also the codimension.");
    }
    if (result.second.second == codim) {
      return result.second.first;
    }
    while (begin->second.second != codim && begin != end) {
      ++begin;
    }
    if (begin->second.second == codim) {
      return begin->second.first;
    }
  }
  throw base::LfException("Physical Entity with name='" + name +
                          "' and codimension=" + std::to_string(codim) +
                          "' not found.");
}

std::string GmshReader::PhysicalEntityNr2Name(size_type number,
                                              dim_t codim) const {
  auto [begin, end] = nr_2_name_.equal_range(number);  // NOLINT
  if (begin == end) {
    throw base::LfException("Physical entity with number " +
                            std::to_string(number) + " not found.");
  }
  auto result = *begin;
  ++begin;
  if (begin == end) {
    if (codim == dim_t(-1) || result.second.second == codim) {
      return result.second.first;
    }
  } else {
    if (codim == dim_t(-1)) {
      throw base::LfException(
          "There are multiple physical entities with the Number " +
          std::to_string(number) + ", please specify also the codimension");
    }
    if (result.second.second == codim) {
      return result.second.first;
    }
    while (begin->second.second != codim && begin != end) {
      ++begin;
    }
    if (begin->second.second == codim) {
      return begin->second.first;
    }
  }
  throw base::LfException(
      "Physical entity with number=" + std::to_string(number) +
      ", codim=" + std::to_string(codim) + " not found.");
}

std::vector<std::pair<size_type, std::string>> GmshReader::PhysicalEntities(
    dim_t codim) const {
  std::vector<std::pair<size_type, std::string>> result;
  for (auto& p : nr_2_name_) {
    if (p.second.second != codim) {
      continue;
    }
    result.emplace_back(p.first, p.second.first);
  }
  return result;
}

std::vector<size_type> GmshReader::PhysicalEntityNr(
    const mesh::Entity& e) const {
  return physical_nrs_->operator()(e);
}

std::variant<GMshFileV2, GMshFileV4> ReadGmshFile(const std::string& filename) {
  // Open file and copy it into memory:
  /////////////////////////////////////////////////////////////////////////////
  std::ifstream in(filename, std::ios_base::in | std::ios_base::binary);
  if (!in) {
    std::string error("Could not open file ");
    error += filename;
    throw lf::base::LfException(error);
  }

  std::string storage;
  in.unsetf(std::ios::skipws);  // no white space skipping
  std::copy(std::istream_iterator<char>(in), std::istream_iterator<char>(),
            std::back_inserter(storage));

  // Parse header to determine if we are dealing with ASCII format or binary
  // format + little or big endian:
  /////////////////////////////////////////////////////////////////////////////
  auto iter = storage.cbegin();
  auto end = storage.cend();

  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;

  // version, is_binary, sizeof(size_t), 1 (as int)
  std::tuple<std::string, bool, int, int> header;
  qi::rule<decltype(iter), std::string()> version_parser;

  version_parser %= +(ascii::alnum | ascii::punct);
  qi::rule<decltype(iter), decltype(header)(), ascii::space_type> header_parser;
  header_parser %= qi::lit("$MeshFormat") >>
                   (qi::hold[(version_parser >> qi::lit('0') >>
                              qi::attr(false) >> qi::int_ >> qi::attr(1))] |
                    (version_parser >> qi::lit('1') >> qi::attr(true) >>
                     qi::int_ >> qi::little_dword)) >>
                   qi::lit("$EndMeshFormat");

  bool succesful;
  succesful = qi::phrase_parse(iter, end, header_parser, ascii::space, header);

  LF_VERIFY_MSG(succesful, "Could not read header of file " << filename);

  if (std::get<0>(header) == "4.1") {
    return ReadGmshFileV4(iter, end, std::get<0>(header), std::get<1>(header),
                          std::get<2>(header), std::get<3>(header), filename);
  } else {
    LF_VERIFY_MSG(false, "GmshFiles with Version "
                             << std::get<0>(header)
                             << " are not yet supported.");
  }
}

}  // namespace lf::io
