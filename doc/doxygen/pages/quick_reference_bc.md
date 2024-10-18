# Quick Reference - Boundary Conditions {#quick_reference_bc}

[TOC]


FixFlaggedSolutionComponents (in both versions) requires a lambda function std::pair<bool,
double> selector(unsigned int dof_idx). It returns wether the dof is to be fixed and if so, the
value it should be fixed to.
We get the boundary flags with
auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2);
The second argument is the codim of entities we are interested in. To fix all dofs on the boundary
(including those associated with edges and cells), give no second argument.
1. One function 𝑔 for the whole boundary
To get the function values of 𝑔 at the boundary nodes, wrap it in a MeshFunction:
auto mf_g = lf::mesh::utils::MeshFunctionGlobal(g);
auto boundary_val =
lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g);
This gives a std::vector<std::pair<bool, scalar_t>>. To get our selector:
auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
return boundary_val[dof_idx];
};
If 𝑔 has diﬀerent definitions on diﬀerent parts of the boundary (e.g., 𝑔 = 0 on Γ0 and 𝑔 = 1 on
Γ1), try to express 𝑔 as a lambda function with an if-else statement.
1. Constant value
The selector becomes
auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
if (bd_flags[dof_idx]) {
return std::make_pair(true, boundary_value);
} else {
return std::make_pair(false, 0.0); // value irrelevant
}
};
1. BC only on part of the boundary
Wee need to create our own bd_flags. Initialize it with default value false:
lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(mesh_p, false);
Then loop over nodes and edges separately and set bd_flags to true for the nodes/edges where
the BC should be applied.
for (const auto& edge : fe_space->Mesh()->Entities(1)) {
if (...) bd_flags(*edge) = true;
}
for (const auto& node : fe_space->Mesh()->Entities(2)) {...}