#ifndef __5ac8f981f27e45d3b9d15fc9d52f7136
#define __5ac8f981f27e45d3b9d15fc9d52f7136

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include "mesh_interface.h"

namespace lf::mesh {

class MeshBuilder {
public:
  /** @copydoc Mesh::size_type */
  using size_type = unsigned int;

  /** @copydoc Geometry::coord_t */
  using coord_t = Eigen::VectorXd;

  using dim_t = unsigned char;


  /**
   * @brief Return the Mesh::DimWorld() of the mesh that will be returned.
   */
  virtual dim_t DimWorld() const = 0;


  /**
   * @brief Return the Mesh::DimMesh() of the mesh that will be returned.
   */
  virtual dim_t DimMesh() const = 0;

  /**
   * @brief Add a point to the mesh.
   * @param coord The coordinate of the point (should have DimWorld() rows)
   * @return The 0-based index of the entity that will be created.
   *         The first call to this method will return 0, the second call 1,
   *         ...
   */
  virtual size_type AddPoint(coord_t coord) = 0;


  /**
   * @brief Add an element (codim=0) to the mesh.
   * @param nodes The 0-based indices of the node that make up this element
   *              (as returned from AddPoint())
   * @param geometry The geometric description of the mesh element.
   * @return The index of the entity that will be created.
   *         The first call to this method will return 0, the second call 1,
   *         and so forth.
   *         
   * @note The node indices passed with the parameter `nodes` must be added
   *       with AddPoint() before calling this method.
   */
  virtual size_type AddElement(const base::ForwardRange<const size_type>& nodes,
                                  std::unique_ptr<geometry::Geometry>&&
                                  geometry) = 0;


  /**
   * @brief Construct a mesh out of the specified nodes and elements.
   * @return The created mesh.
   * 
   * @note After calling Build() you are not allowed to call AddPoint or 
   *       AddElement() anymore.
   */
  virtual std::unique_ptr<Mesh> Build() = 0;

  /// @brief Virtual destructor.
  virtual ~MeshBuilder() = default;
};

} // namespace lf::mesh

#endif  // __5ac8f981f27e45d3b9d15fc9d52f7136
