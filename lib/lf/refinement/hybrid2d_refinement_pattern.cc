/**
 * @file refinement.cc
 * @brief implementation of splitting of reference entities according to
 * refinement patterns
 */

#include "hybrid2d_refinement_pattern.h"

namespace lf::refinement {

std::ostream &operator<<(std::ostream &o, const RefPat &refpat) {
  switch (refpat) {
  rp_nil : {
    o << "rp_NIL";
    break;
  }
  rp_copy : {
    o << "rp_COPY";
    break;
  }
  rp_split : {
    o << "rp_SPLIT";
    break;
  }
  rp_bisect : {
    o << "rp_BISECT";
    break;
  }
  rp_trisect : {
    o << "rp_TRISECT";
    break;
  }
  rp_trisect_left : {
    o << "rp_TRISECT_LEFT";
    break;
  }
  rp_quadsect : {
    o << "rp_QUADSECT";
    break;
  }
  rp_threeedge : {
    o << "rp_THREEEDGE";
    break;
  }
  rp_regular : {
    o << "rp_REGULAR";
    break;
  }
  rp_barycentric : {
    o << "rp_BARYCENTRIC";
    break;
  }
  case rp_nil:
    o << "rp_nil";
    break;
    case rp_copy:
      o << "rp_copy";
      break;
    case rp_split:
      o << "rp_split";
      break;
    case rp_bisect:
      o << "rp_bisect";
      break;
    case rp_trisect:
      o << "rp_trisect";
      break;
    case rp_trisect_left:
      o << "rp_trisect_left";
      break;
    case rp_quadsect:
      o << "rp_quadsect";
      break;
    case rp_threeedge:
      o << "rp_threeedge";
      break;
    case rp_regular:
      o << "rp_regular";
      break;
    case rp_barycentric:
      o << "rp_barycentric";
      break;
  }
  return o;
}

// Implementation for RefinemeentPattern

lf::base::size_type Hybrid2DRefinementPattern::NumChildren(
    lf::base::dim_t codim) const {
  LF_VERIFY_MSG(codim <= ref_el_.Dimension(),
                "Illegal codim = " << codim << " for " << ref_el_.ToString());
  // Depending on the type of cell do something different
  switch (ref_el_) {
    case lf::base::RefEl::kPoint(): {
      switch (ref_pat_) {
        case RefPat::rp_nil: {
          return 0;
        }
        case RefPat::rp_copy: {
          return 1;
        }
        default: {
          LF_VERIFY_MSG(false, "Illegal refinement pattern for point");
        }
      }  // end swtich ref_pat
      break;
    }  // end case of a simple point
    case lf::base::RefEl::kSegment(): {
      switch (codim) {
        case 0: {
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              return 0;
            }
            case RefPat::rp_copy: {
              return 1;
            }
            case RefPat::rp_split: {
              return 2;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for edge!");
              break;
            }
          }  // end switch ref_pat
          break;
        }
        case 1: {
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              return 0;
            }
            case RefPat::rp_copy: {
              return 0;
            }
            case RefPat::rp_split: {
              return 1;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for edge!");
              break;
            }
          }  // end switch ref_pat
          break;
        }
        default:
          LF_VERIFY_MSG(false, "invalid codim");
      }  // end switch codim
      break;
    }  // end case of an edge
    case lf::base::RefEl::kTria(): {
      switch (codim) {
        case 0: {
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              return 0;
            }
            case RefPat::rp_copy: {
              return 1;
            }
            case RefPat::rp_bisect: {
              return 2;
            }
            case RefPat::rp_trisect: {
              return 3;
            }
            case RefPat::rp_trisect_left: {
              return 3;
            }
            case RefPat::rp_quadsect: {
              return 4;
            }
            case RefPat::rp_regular: {
              return 4;
            }
            case RefPat::rp_barycentric: {
              return 6;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for triangle!");
              break;
            }
          }  // end switch ref_pat
          break;
          default:
            LF_VERIFY_MSG(false, "invalid codim");
        }  // end case codim = 0
        case 1: {
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              return 0;
            }
            case RefPat::rp_copy: {
              return 0;
            }
            case RefPat::rp_bisect: {
              return 1;
            }
            case RefPat::rp_trisect: {
              return 2;
            }
            case RefPat::rp_trisect_left: {
              return 2;
            }
            case RefPat::rp_quadsect: {
              return 3;
            }
            case RefPat::rp_regular: {
              return 3;
            }
            case RefPat::rp_barycentric: {
              return 6;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << ref_pat_
                                       << "not implemented for triangle!");
              break;
            }
          }  // end switch ref_pat
          break;
        }  // end case codim = 1
        case 2: {
          if (ref_pat_ == RefPat::rp_barycentric) {
            return 1;
          }
          return 0;
        }
      }  // end switch codim
      break;
    }  // end case of a triangle
    case lf::base::RefEl::kQuad(): {
      switch (codim) {
        case 0: {
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              return 0;
            }
            case RefPat::rp_copy: {
              return 1;
            }
            case RefPat::rp_trisect: {
              return 3;
            }
            case RefPat::rp_quadsect: {
              return 4;
            }
            case RefPat::rp_bisect:
            case RefPat::rp_split: {
              return 2;
            }
            case RefPat::rp_threeedge: {
              return 4;
            }
            case RefPat::rp_barycentric:
            case RefPat::rp_regular: {
              return 4;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << ref_pat_
                                       << "not implemented for quadrilateral!");
              break;
            }
          }  // end switch ref_pat
          break;
        }
        case 1: {
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              return 0;
            }
            case RefPat::rp_copy: {
              return 0;
            }
            case RefPat::rp_trisect: {
              return 2;
            }
            case RefPat::rp_quadsect: {
              return 3;
            }
            case RefPat::rp_bisect:
            case RefPat::rp_split: {
              return 1;
            }
            case RefPat::rp_threeedge: {
              return 3;
            }
            case RefPat::rp_barycentric:
            case RefPat::rp_regular: {
              return 4;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for quadrilateral!");
              break;
            }
          }  // end switch ref_pat
          break;
        }
        case 2: {
          if ((ref_pat_ == RefPat::rp_barycentric) ||
              (ref_pat_ == RefPat::rp_regular)) {
            return 1;
          }
          return 0;
          break;
        }
        default:
          LF_VERIFY_MSG(false, "invalid codim");
      }  // end switch codim
      break;
    }  // end case of a quadrilateral
  }    // end switch cell type
  return 0;
}

std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
Hybrid2DRefinementPattern::ChildPolygons(lf::base::dim_t codim) const {
  LF_VERIFY_MSG(codim <= ref_el_.Dimension(),
                "Illegal codim = " << codim << " for " << ref_el_.ToString());
  // Local variable for accumulating information about interior entities
  // created during refinement
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> child_poly{};

  // Lattice point coordinates (lattice_const_ should be a multiple of 6)
  const unsigned int lt_half = lattice_const_ / 2;
  const unsigned int lt_third = lattice_const_ / 3;
  const unsigned int lt_one = lattice_const_;
  // Depending on the type of cell do something different
  switch (ref_el_) {
    case lf::base::RefEl::kPoint(): {
      if (ref_pat_ != RefPat::rp_nil) {
        child_poly.emplace_back(0, 1);
      }
      break;
    }                                    // end case of a simple point
    case lf::base::RefEl::kSegment(): {  // case of an edge
      switch (codim) {
        case 0: {  // child entities are segments as well
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              break;
            }
            case RefPat::rp_copy: {  // copy of segment, 1 child
              Eigen::Matrix<int, 1, 2> seg_ref_coord;
              seg_ref_coord << 0, lt_one;
              child_poly.emplace_back(seg_ref_coord);
              break;
            }
            case RefPat::rp_split: {  // splitting of segment, 2 children
              Eigen::Matrix<int, 1, 2> part_ref_coord;
              part_ref_coord << 0, lt_half;
              child_poly.emplace_back(part_ref_coord);
              part_ref_coord << lt_half, lt_one;
              child_poly.emplace_back(part_ref_coord);
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for edge!");
              break;
            }
          }  // end switch ref_pat
          break;
        }          // case codim = 0
        case 1: {  // child entities are points
          if (ref_pat_ == RefPat::rp_split) {
            // child is the midpoint
            Eigen::Matrix<int, 1, 1> point_ref_coord;
            point_ref_coord << lt_half;
            child_poly.emplace_back(point_ref_coord);
          }
          break;
        }  // end case codim = 1
        default:
          LF_VERIFY_MSG(false, "invalid codim");
      }  // end switch codim
      break;
    }  // end case of an edge
    case lf::base::RefEl::kTria(): {
      // **********************************************************************
      // Triangle
      // **********************************************************************
      // Lattice coordinates of special points in the triangle
      // Here we use knowledge about the conventions underlying node
      // and edge numbering in a triangle as defined in RefEl.
      // This information could also be retrieved from the reference element.
      Eigen::MatrixXi lt_node_coords(2, 3);
      lt_node_coords << 0, lt_one, 0, 0, 0, lt_one;
      Eigen::MatrixXi lt_midpoint_coords(2, 3);
      lt_midpoint_coords << lt_half, lt_half, 0, 0, lt_half, lt_half;

      // Remap local indices according to anchor values
      const unsigned int mod_0 = (0 + anchor_) % 3;
      const unsigned int mod_1 = (1 + anchor_) % 3;
      const unsigned int mod_2 = (2 + anchor_) % 3;

      switch (codim) {
        case 0: {
          // For temporary storage of triangular lattice polygon
          Eigen::MatrixXi child_coords(2, 3);

          switch (ref_pat_) {
            case RefPat::rp_nil: {
              break;
            }
            case RefPat::rp_copy: {
              child_coords = lt_node_coords;
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_bisect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for bisection refinement of triangle");
              // Splitting a triangle in two by bisecting anchor edge
              child_coords.col(0) = lt_node_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_1);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_trisect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for trisection refinement of triangle");
              // Bisect through anchor edge first and then bisect through
              // edge with the next larger index (mod 3); creates three
              // child triangles.
              child_coords.col(0) = lt_node_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_1);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_2);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_trisect_left: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for trisection refinement of triangle");

              // Bisect through anchor edge first and then bisect through
              // edge with the next smaller index (mod 3); creates three
              // child triangles.
              child_coords.col(0) = lt_node_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_1);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_2);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_quadsect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for quadsection refinement of triangle");
              // Bisect through the anchor edge first and then
              // through the two remaining edges; creates four child
              // triangles; every edge is split.
              child_coords.col(0) = lt_node_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_1);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_2);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_2);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_coords.col(2) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_regular: {
              // Split triangle into four small  congruent triangles
              child_coords.col(0) = lt_node_coords.col(0);
              child_coords.col(1) = lt_midpoint_coords.col(0);
              child_coords.col(2) = lt_midpoint_coords.col(2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(1);
              child_coords.col(1) = lt_midpoint_coords.col(0);
              child_coords.col(2) = lt_midpoint_coords.col(1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(2);
              child_coords.col(1) = lt_midpoint_coords.col(2);
              child_coords.col(2) = lt_midpoint_coords.col(1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(0);
              child_coords.col(1) = lt_midpoint_coords.col(1);
              child_coords.col(2) = lt_midpoint_coords.col(2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_barycentric: {
              //  Split triangle into 5 smaller triangles by connecting
              // the center of gravity with the vertices and the midpoints
              // of the edges.

              // Obtain coordinates of center of gravity
              Eigen::Vector2i lt_baryc_coords =
                  Eigen::Vector2i({lt_third, lt_third});

              child_coords.col(0) = lt_node_coords.col(0);
              child_coords.col(1) = lt_midpoint_coords.col(0);
              child_coords.col(2) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(1);
              child_coords.col(1) = lt_midpoint_coords.col(0);
              child_coords.col(2) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(1);
              child_coords.col(1) = lt_midpoint_coords.col(1);
              child_coords.col(2) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(2);
              child_coords.col(1) = lt_midpoint_coords.col(1);
              child_coords.col(2) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(2);
              child_coords.col(1) = lt_midpoint_coords.col(2);
              child_coords.col(2) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(0);
              child_coords.col(1) = lt_midpoint_coords.col(2);
              child_coords.col(2) = lt_baryc_coords;
              child_poly.push_back(child_coords);
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for triangle!");
              break;
            }
          }  // end switch ref_pat
          break;
        }          // end case codim = 0
        case 1: {  // codim == 1
          // Return the location of interior edges of the triangle
          // This location is defined by two lattice points given by the
          // rows of the following matrix
          Eigen::MatrixXi child_coords(2, 2);

          // Different layout of interior edges for different refinement
          // patterns
          switch (ref_pat_) {
            case RefPat::rp_nil:
            case RefPat::rp_copy: {
              // No interior edges in these cases
              break;
            }
            case RefPat::rp_bisect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for bisection refinement of triangle");
              // Splitting a triangle in two by bisecting anchor edge
              // One interior edge is created
              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_trisect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for trisection refinement of triangle");
              // Bisect through anchor edge first and then bisect through
              // edge with the next larger index (mod 3); creates two
              // interior edges
              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_trisect_left: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for trisection refinement of triangle");

              // Bisect through anchor edge first and then bisect through
              // edge with the next smaller index (mod 3); creates two
              // interior edges
              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_quadsect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for quadsection refinement of triangle");
              // Bisect through the anchor edge first and then
              // through the two remaining edges; creates three
              // interior edges
              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_regular: {
              // Split triangle into four small  congruent triangles
              // Introduces three interior edges
              child_coords.col(0) = lt_midpoint_coords.col(0);
              child_coords.col(1) = lt_midpoint_coords.col(2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(0);
              child_coords.col(1) = lt_midpoint_coords.col(1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(2);
              child_coords.col(1) = lt_midpoint_coords.col(1);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_barycentric: {
              //  Split triangle into 5 smaller triangles by connecting
              // the center of gravity with the vertices and the midpoints
              // of the edges. Six interior edges will be created

              // Obtain coordinates of center of gravity
              Eigen::Vector2i lt_baryc_coords =
                  Eigen::Vector2i({lt_third, lt_third});

              child_coords.col(0) = lt_node_coords.col(0);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(1);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(2);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(0);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(1);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(2);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for triangle!");
              break;
            }
          }  // end switch ref_pat
          break;
        }  // end codim = 1
        case 2: {
          // Return interior points
          if (ref_pat_ == RefPat::rp_barycentric) {
            // Only in the case of barycentric refinemnt insert new point
            // in the center of gravitfy.
            child_poly.emplace_back(Eigen::Vector2i({lt_third, lt_third}));
          }
          break;
        }  // end codim == 2
        default:
          LF_VERIFY_MSG(false, "invalid codim");
      }  // end switch codim
      break;
    }  // end case of a triangle
    case lf::base::RefEl::kQuad(): {
      // Lattice coordinates of special points in a quadrilateral
      // Here we use knowledge about the conventions underlying node
      // and edge numbering in a quadrilateral as defined in RefEl.
      // This information could also be retrieved from the reference element.
      Eigen::MatrixXi lt_node_coords(2, 4);
      lt_node_coords << 0, lt_one, lt_one, 0, 0, 0, lt_one, lt_one;
      Eigen::MatrixXi lt_midpoint_coords(2, 4);
      lt_midpoint_coords << lt_half, lt_one, lt_half, 0, 0, lt_half, lt_one,
          lt_half;

      // Remap local indices according to anchor values
      const unsigned int mod_0 = (0 + anchor_) % 4;
      const unsigned int mod_1 = (1 + anchor_) % 4;
      const unsigned int mod_2 = (2 + anchor_) % 4;
      const unsigned int mod_3 = (3 + anchor_) % 4;

      switch (codim) {
        case 0: {
          // Return sub-cell location inside quadrilateral
          // For temporary storage of triangular lattice polygon
          Eigen::MatrixXi tria_child_coords(2, 3);
          // For temporary storage of 4-node lattice polygon
          Eigen::MatrixXi quad_child_coords(2, 4);
          switch (ref_pat_) {
            case RefPat::rp_nil: {
              break;
            }
            case RefPat::rp_copy: {
              child_poly.push_back(lt_node_coords);
              break;
            }
            case RefPat::rp_trisect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for trisection refinement of quad");

              // Partition a quad into three triangle, the anchor edge
              // being split in the process
              tria_child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(1) = lt_node_coords.col(mod_2);
              tria_child_coords.col(2) = lt_node_coords.col(mod_3);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(1) = lt_node_coords.col(mod_0);
              tria_child_coords.col(2) = lt_node_coords.col(mod_3);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(1) = lt_node_coords.col(mod_1);
              tria_child_coords.col(2) = lt_node_coords.col(mod_2);
              child_poly.push_back(tria_child_coords);
              break;
            }
            case RefPat::rp_quadsect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for quadsection refinement of triangle");
              // Partition a quad into four triangle, thus
              // splitting two edges. The one with the smaller sub index is the
              // anchor edge
              tria_child_coords.col(0) = lt_node_coords.col(mod_0);
              tria_child_coords.col(1) = lt_node_coords.col(mod_3);
              tria_child_coords.col(2) = lt_midpoint_coords.col(mod_0);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_node_coords.col(mod_1);
              tria_child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              tria_child_coords.col(2) = lt_midpoint_coords.col(mod_0);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_node_coords.col(mod_2);
              tria_child_coords.col(1) = lt_node_coords.col(mod_3);
              tria_child_coords.col(2) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              tria_child_coords.col(2) = lt_node_coords.col(mod_3);
              child_poly.push_back(tria_child_coords);
              break;
            }
            case RefPat::rp_bisect:
            case RefPat::rp_split: {
              LF_VERIFY_MSG(anchor_set_,
                            "Anchor must be set for splitting of quad");

              // Cut a quadrilateral into two
              quad_child_coords.col(0) = lt_node_coords.col(mod_0);
              quad_child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              quad_child_coords.col(2) = lt_midpoint_coords.col(mod_2);
              quad_child_coords.col(3) = lt_node_coords.col(mod_3);
              child_poly.push_back(quad_child_coords);

              quad_child_coords.col(0) = lt_node_coords.col(mod_1);
              quad_child_coords.col(1) = lt_node_coords.col(mod_2);
              quad_child_coords.col(2) = lt_midpoint_coords.col(mod_2);
              quad_child_coords.col(3) = lt_midpoint_coords.col(mod_0);
              child_poly.push_back(quad_child_coords);
              break;
            }
            case RefPat::rp_threeedge: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for three edge refinement of a quad");

              // A quadrilateral with three split edges is decomposed into one
              // child quadrilateral and three triangles
              quad_child_coords.col(0) = lt_node_coords.col(mod_2);
              quad_child_coords.col(1) = lt_node_coords.col(mod_3);
              quad_child_coords.col(2) = lt_midpoint_coords.col(mod_3);
              quad_child_coords.col(3) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(quad_child_coords);

              tria_child_coords.col(0) = lt_node_coords.col(mod_0);
              tria_child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(2) = lt_midpoint_coords.col(mod_3);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_node_coords.col(mod_1);
              tria_child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(2) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(tria_child_coords);

              tria_child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              tria_child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              tria_child_coords.col(2) = lt_midpoint_coords.col(mod_3);
              child_poly.push_back(tria_child_coords);
              break;
            }
            case RefPat::rp_barycentric:
            case RefPat::rp_regular: {
              // Fully symmetric splitting into four quadrilaterals
              // Obtain coordinates of center of gravity
              Eigen::Vector2i lt_baryc_coords =
                  Eigen::Vector2i({lt_half, lt_half});

              quad_child_coords.col(0) = lt_node_coords.col(0);
              quad_child_coords.col(1) = lt_midpoint_coords.col(0);
              quad_child_coords.col(2) = lt_baryc_coords;
              quad_child_coords.col(3) = lt_midpoint_coords.col(3);
              child_poly.push_back(quad_child_coords);

              quad_child_coords.col(0) = lt_node_coords.col(1);
              quad_child_coords.col(1) = lt_midpoint_coords.col(1);
              quad_child_coords.col(2) = lt_baryc_coords;
              quad_child_coords.col(3) = lt_midpoint_coords.col(0);
              child_poly.push_back(quad_child_coords);

              quad_child_coords.col(0) = lt_node_coords.col(2);
              quad_child_coords.col(1) = lt_midpoint_coords.col(1);
              quad_child_coords.col(2) = lt_baryc_coords;
              quad_child_coords.col(3) = lt_midpoint_coords.col(2);
              child_poly.push_back(quad_child_coords);

              quad_child_coords.col(0) = lt_node_coords.col(3);
              quad_child_coords.col(1) = lt_midpoint_coords.col(2);
              quad_child_coords.col(2) = lt_baryc_coords;
              quad_child_coords.col(3) = lt_midpoint_coords.col(3);
              child_poly.push_back(quad_child_coords);
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for quadrilateral!");
              break;
            }
          }  // end switch ref_pat
          break;
        }  // end codim == 0
        case 1: {
          // Return information about edges created inside a quadrilateral
          Eigen::MatrixXi child_coords(2, 2);
          switch (ref_pat_) {
            case RefPat::rp_nil:
            case RefPat::rp_copy: {
              // No internal edges in this case
              break;
            }
            case RefPat::rp_trisect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for trisection refinement of quad");

              // Partition a quad into three triangle, the anchor edge
              // being split in the process. Two interior edges are created
              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_node_coords.col(mod_2);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_node_coords.col(mod_3);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_quadsect: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for quadsection refinement of triangle");
              // Partition a quad into four triangle, thus
              // splitting two edges. This spawns three internal edges
              child_coords.col(0) = lt_node_coords.col(mod_3);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_1);
              child_coords.col(1) = lt_midpoint_coords.col(mod_0);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_node_coords.col(mod_3);
              child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_bisect:
            case RefPat::rp_split: {
              LF_VERIFY_MSG(anchor_set_,
                            "Anchor must be set for splitting of quad");

              // Cut a quadrilateral into two.  One interior edge arises
              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_2);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_threeedge: {
              LF_VERIFY_MSG(
                  anchor_set_,
                  "Anchor must be set for three edge refinement of a quad");
              // Split quadrilateral into one half and three triangles, creating
              // three interior edges in the process
              child_coords.col(0) = lt_midpoint_coords.col(mod_3);
              child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_3);
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(mod_0);
              child_coords.col(1) = lt_midpoint_coords.col(mod_1);
              child_poly.push_back(child_coords);
              break;
            }
            case RefPat::rp_barycentric:
            case RefPat::rp_regular: {
              // Fully symmetric splitting into four quadrilaterals
              // Four interior edges arise
              // Obtain coordinates of center of gravity
              Eigen::Vector2i lt_baryc_coords =
                  Eigen::Vector2i({lt_half, lt_half});

              child_coords.col(0) = lt_midpoint_coords.col(0);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(1);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(2);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);

              child_coords.col(0) = lt_midpoint_coords.col(3);
              child_coords.col(1) = lt_baryc_coords;
              child_poly.push_back(child_coords);
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "refinement pattern "
                                       << (int)ref_pat_
                                       << "not implemented for quadrilateral!");
              break;
            }
          }  // end switch ref_pat
          break;
        }
        case 2: {
          // Return location of new interior points
          // Those arise only in the case of regular/barycentric refinement
          if ((ref_pat_ == RefPat::rp_regular) ||
              (ref_pat_ == RefPat::rp_barycentric)) {
            child_poly.emplace_back(Eigen::Vector2i({lt_half, lt_half}));
          }
          break;
        }  // end codim == 2
        default:
          LF_VERIFY_MSG(false, "invalid codim");
      }  // end switch codim
      break;
    }  // end case of a quadrilateral
  }    // end switch cell type
  return (child_poly);
}

}  // namespace lf::refinement
