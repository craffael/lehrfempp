% MATLAB script for ploting testing refinement hierarchy

display('#### Regularly refined test mesh ####');
opts = []; opts.numbers = true; plot_lf_mesh(@fine_mesh,opts);
title('Regularly refined test mesh');
print -dpdf -bestfit 'fine_mesh.pdf'
opts = []; opts.parents = @fine_mesh_pi; plot_lf_mesh(@fine_mesh,opts);
title('Regularly refined test mesh: parent info');
print -dpdf -bestfit 'fine_mesh_pi.pdf'

display('#### Barycentric refinement of test mesh ####');
opts = []; opts.numbers = true; plot_lf_mesh(@barycentric_ref,opts);
title('Barycentric refinement of test mesh');
print -dpdf -bestfit 'barycentric_ref.pdf'
opts = []; opts.parents = @barycentric_ref_pi; plot_lf_mesh(@barycentric_ref,opts);
title('Barycentric refinement of test mesh: parent info');
print -dpdf -bestfit 'barycentric_ref_pi.pdf'

display('#### Global bisection refinement of test mesh ####');
opts = []; opts.numbers = true; plot_lf_mesh(@allref,opts);
title('Bisection refinement of test mesh');
print -dpdf -bestfit 'allref.pdf'
opts = []; opts.parents = @allref_pi; plot_lf_mesh(@allref,opts);
title('Bisection refinement of test mesh: parent info');
print -dpdf -bestfit 'allref_pi.pdf'

display('#### Hierarchy of locally refined meshes ######');

display('Level 0');
opts = []; opts.numbers = true; plot_lf_mesh(@multiref_L0,opts);
title('Test mesh (= level 0)');
print -dpdf -bestfit 'multiref_L0.pdf'

display('Level 1');
opts = []; ; plot_lf_mesh(@multiref_L1,opts);
title(sprintf('Multirev level %i',1));
print -dpdf -bestfit 'multiref_L1.pdf'
opts = []; opts.parents = @multiref_L1_pi; plot_lf_mesh(@multiref_L1,opts);
title(sprintf('Multirev parent info %i',1));
print -dpdf -bestfit 'multiref_L1_pi.pdf'

display('Level 2');
opts = []; ; plot_lf_mesh(@multiref_L2,opts);
title(sprintf('Multirev level %i',2));
print -dpdf -bestfit 'multiref_L2.pdf'
opts = []; opts.parents = @multiref_L2_pi; plot_lf_mesh(@multiref_L2,opts);
title(sprintf('Multirev parent info %i',2));
print -dpdf -bestfit 'multiref_L2_pi.pdf'

display('Level 3');
opts = []; ; plot_lf_mesh(@multiref_L3,opts);
title(sprintf('Multirev level %i',3));
print -dpdf -bestfit 'multiref_L3.pdf'
opts = []; opts.parents = @multiref_L3_pi; plot_lf_mesh(@multiref_L3,opts);
title(sprintf('Multirev parent info %i',3));
print -dpdf -bestfit 'multiref_L3_pi.pdf'

display('Level 4');
opts = []; ; plot_lf_mesh(@multiref_L4,opts);
title(sprintf('Multirev level %i',4));
print -dpdf -bestfit 'multiref_L4.pdf'
opts = []; opts.parents = @multiref_L4_pi; plot_lf_mesh(@multiref_L4,opts);
title(sprintf('Multirev parent info %i',4));
print -dpdf -bestfit 'multiref_L4_pi.pdf'

display('#### Mixed local/global refinement ######');

display('Level 0');
opts = []; ; plot_lf_mesh(@mixedref_0_L0,opts);
title(sprintf('Multirev level %i',0));
print -dpdf -bestfit 'mixedref_0_L0.pdf'

display('Level 1');
opts = []; ; plot_lf_mesh(@mixedref_0_L1,opts);
title(sprintf('Multirev level %i',1));
print -dpdf -bestfit 'mixedref_0_L1.pdf'
opts = []; opts.parents = @mixedref_0_L1_pi; plot_lf_mesh(@mixedref_0_L1,opts);
title(sprintf('Multirev parent info %i',1));

display('Level 2');
opts = []; ; plot_lf_mesh(@mixedref_0_L2,opts);
title(sprintf('Multirev level %i',2));
print -dpdf -bestfit 'mixedref_0_L2.pdf'
opts = []; opts.parents = @mixedref_0_L2_pi; plot_lf_mesh(@mixedref_0_L2,opts);
title(sprintf('Multirev parent info %i',2));

display('Level 3');
opts = []; ; plot_lf_mesh(@mixedref_0_L3,opts);
title(sprintf('Multirev level %i',3));
print -dpdf -bestfit 'mixedref_0_L3.pdf'
opts = []; opts.parents = @mixedref_0_L3_pi; plot_lf_mesh(@mixedref_0_L3,opts);
title(sprintf('Multirev parent info %i',3));

display('Level 4');
opts = []; ; plot_lf_mesh(@mixedref_0_L4,opts);
title(sprintf('Multirev level %i',4));
print -dpdf -bestfit 'mixedref_0_L4.pdf'
opts = []; opts.parents = @mixedref_0_L4_pi; plot_lf_mesh(@mixedref_0_L4,opts);
title(sprintf('Multirev parent info %i',4));
