#ifndef __e98a803fac5b430a8ff634ceb2f809a1
#define __e98a803fac5b430a8ff634ceb2f809a1

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements mesh::MeshBuilder interface and can be used to construct
 *        a hybrid mesh with `dimMesh=2`.
 */
class MeshBuilder : public mesh::MeshBuilder {
  
public:

  /**
   * @brief Construct a new builder that can be used to construct a new hybrid2d
   *        mesh.
   * @param dim_world The dimension of the euclidean space in which the
   *                  mesh is embedded.
   */
  MeshBuilder(dim_t dim_world) : dim_world_(dim_world), built_(false) {}

  dim_t DimWorld() const override {return dim_world_; }
  dim_t DimMesh() const override {return 2; }

  size_type AddPoint(coord_t coord) override;
  size_type AddElement(const base::ForwardRange<const size_type>& nodes,
    std::unique_ptr<geometry::Geometry>&& geometry) override;
  std::unique_ptr<mesh::Mesh> Build() override;
private:
  dim_t dim_world_;
  bool built_;
  std::vector<Eigen::VectorXd> nodes_;
  std::vector<std::tuple<std::vector<size_type>, std::unique_ptr<geometry::Geometry>>> elements_;
};

}



#endif // __e98a803fac5b430a8ff634ceb2f809a1
