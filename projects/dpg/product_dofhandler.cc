#include "product_dofhandler.h"

//----------------------------------------------
// Implementation of ProductUniformFEDofHandler
//-----------------------------------------------
namespace projects::dpg {

ProductUniformFEDofHandler::ProductUniformFEDofHandler(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    std::vector<dof_map_t> dofmaps)
    : mesh_p_(std::move(mesh_p)), num_components_(dofmaps.size()) {
  // check, that there local shape function layouts are actually specified
  LF_ASSERT_MSG(num_components_ != 0, "dofmaps.size() = 0");

  component_dof_handlers_.resize(num_components_);
  offsets_.resize(num_components_ + 1);
  glb_idx_t offset = 0;

  // initialize the component dof handlers from the specified local shape
  // function layout for each of the components. Keep track of the offsets.
  for (size_type component = 0; component < num_components_; ++component) {
    component_dof_handlers_[component] =
        std::make_shared<lf::assemble::UniformFEDofHandler>(mesh_p_,
                                                            dofmaps[component]);
    offsets_[component] = offset;
    offset += component_dof_handlers_[component]->NumDofs();
  }
  // offsets_[num_components_] =  NoDofs()
  offsets_[num_components_] = offset;

  initIndexArrays();
}

void ProductUniformFEDofHandler::initIndexArrays() {
  // construct the index arrays based on the offsets, the component dofhandlers
  // and the specified access pattern [codim][entity_idx][local_idx].
  // iterate over all codimensions
  for (int codim = 0; codim <= 2; codim++) {
    dofs_[codim].resize(mesh_p_->NumEntities(codim));
    internal_dofs_[codim].resize(mesh_p_->NumEntities(codim));
    // iterate over all entities of a given codimension
    for (glb_idx_t entity_idx = 0; entity_idx < mesh_p_->NumEntities(codim);
         entity_idx++) {
      // retrive entity of a certain index
      const lf::mesh::Entity* entity_p{
          mesh_p_->EntityByIndex(codim, entity_idx)};
      // set up the index mapping using the component dofhandlers and
      // the specified offsets:
      for (size_type component = 0; component < num_components_; component++) {
        // set up local dof indieces mapping
        for (glb_idx_t gdof_idx :
             component_dof_handlers_[component]->GlobalDofIndices(*entity_p)) {
          dofs_[codim][entity_idx].push_back(gdof_idx + offsets_[component]);
        }
        // set up internal dof indices mapping:
        for (glb_idx_t gdof_idx :
             component_dof_handlers_[component]->InteriorGlobalDofIndices(
                 *entity_p)) {
          internal_dofs_[codim][entity_idx].push_back(gdof_idx +
                                                      offsets_[component]);
        }
      }
    }
  }
}

size_type ProductUniformFEDofHandler::NumDofs() const {
  return offsets_[num_components_];
}

size_type ProductUniformFEDofHandler::NumDofs(size_type component) const {
  return component_dof_handlers_[component]->NumDofs();
}

size_type ProductUniformFEDofHandler::NumLocalDofs(
    const lf::mesh::Entity& entity) const {
  const dim_t codim = 2 - entity.RefEl().Dimension();
  const glb_idx_t entity_idx = mesh_p_->Index(entity);
  size_type num_dofs = dofs_[codim][entity_idx].size();
  return num_dofs;
}

size_type ProductUniformFEDofHandler::NumLocalDofs(
    const lf::mesh::Entity& entity, size_type component) const {
  return component_dof_handlers_[component]->NumLocalDofs(entity);
}

size_type ProductUniformFEDofHandler::NumInteriorDofs(
    const lf::mesh::Entity& entity) const {
  const dim_t codim = 2 - entity.RefEl().Dimension();
  const glb_idx_t entity_idx = mesh_p_->Index(entity);
  size_type num_dofs = internal_dofs_[codim][entity_idx].size();
  return num_dofs;
}

size_type ProductUniformFEDofHandler::NumInteriorDofs(
    const lf::mesh::Entity& entity, size_type component) const {
  return component_dof_handlers_[component]->NumInteriorDofs(entity);
}

nonstd::span<const gdof_idx_t> ProductUniformFEDofHandler::GlobalDofIndices(
    const lf::mesh::Entity& entity) const {
  const dim_t codim = 2 - entity.RefEl().Dimension();
  const glb_idx_t entity_idx = mesh_p_->Index(entity);

  const gdof_idx_t* begin = dofs_[codim][entity_idx].data();
  const gdof_idx_t* end = begin + dofs_[codim][entity_idx].size();

  return {begin, end};
}

nonstd::span<const gdof_idx_t>
ProductUniformFEDofHandler::InteriorGlobalDofIndices(
    const lf::mesh::Entity& entity) const {
  const dim_t codim = 2 - entity.RefEl().Dimension();
  const glb_idx_t entity_idx = mesh_p_->Index(entity);

  const gdof_idx_t* begin = internal_dofs_[codim][entity_idx].data();
  const gdof_idx_t* end = begin + internal_dofs_[codim][entity_idx].size();
  return {begin, end};
}

std::shared_ptr<const lf::mesh::Mesh> ProductUniformFEDofHandler::Mesh() const {
  return mesh_p_;
}

const lf::mesh::Entity& ProductUniformFEDofHandler::Entity(
    gdof_idx_t dofnum) const {
  size_type component = Component(dofnum);
  glb_idx_t offset = offsets_[component];
  return component_dof_handlers_[component]->Entity(dofnum - offset);
}

ldof_idx_t ProductUniformFEDofHandler::LocalStartIndex(
    const lf::mesh::Entity& entity, size_type component) const {
  size_type local_offset = 0;
  for (size_type comp = 0; comp < component; comp++) {
    local_offset += NumLocalDofs(entity, comp);
  }
  return local_offset;
}

size_type ProductUniformFEDofHandler::Component(glb_idx_t dofnum) const {
  LF_ASSERT_MSG(dofnum < NumDofs(), "invalid global index (out of bounds)");
  // iterate through compont offsets:
  for (size_type component = 0; component < num_components_; component++) {
    if (offsets_[component + 1] > dofnum) {
      return component;
    }
  }
  LF_ASSERT_MSG(false, "invalid dof number ");
  return 0;
}

}  // namespace projects::dpg
