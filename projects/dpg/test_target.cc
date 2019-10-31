
#include <lf/mesh/test_utils/test_meshes.h>
#include <iostream>

#include "discontinuous_fe_constant.h"
#include "discontinuous_scalar_reference_finite_element.h"
#include "lagr_fe_cubic.h"
#include "lagr_fe_quadratic.h"
#include "trace_scalar_reference_finite_element.h"

#include "product_dofhandler.h"
#include "product_fe_space.h"
#include "product_fe_space_factory.h"

#include "sub_element_matrix_provider.h"
#include "sub_element_vector_provider.h"

#include "product_element_matrix_provider.h"
#include "product_element_matrix_provider_factory.h"
#include "product_element_vector_provider.h"
#include "product_element_vector_provider_factory.h"

#include "loc_comp_dpg.h"

#include "dpg_tools.h"

namespace projects::dpg {

int main() {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  FeLagrangeO3Tria<double> tfe{};
  /*
    ProductUniformFESpaceFactory<double> fac(mesh_p);
    auto u = fac.AddH1Component(2);
    auto fe_space_trial = fac.Build();


    auto mf_one = lf::uscalfe::MeshFunctionConstant(1.0);
    ProductElementMatrixProviderFactory fact(fe_space_trial,fe_space_trial);
    fact.AddDiffusionElementMatrixProvider(u,u,mf_one);

    auto provider = fact.Build();

    for(const auto* const cell: mesh_p->Entities(0)){
      std::cout << provider->Eval(*cell) << std::endl;
    }

    std::cout << "sucessfull exectution of test target \n";
    return 0;
    */
  return 0;
}
}  // namespace projects::dpg
