% MATLAB driver script for plotting the output of LehreFEM++ 
% demo code lf.examples.mesh.test_mesh_output_demo 

display('Plotting of meshes output by lf.examples.mesh.test_mesh_output_demo');
display('(Utility function plot_lf_mesh() must be available)');

opts.numbers = true;

display('Test mesh 0');
plot_lf_mesh(@test_mesh_0,opts);
print -dpng 'test_mesh_0.png';

display('Test mesh 1');
plot_lf_mesh(@test_mesh_1,opts);
print -dpng 'test_mesh_1.png';

display('Test mesh 2');
plot_lf_mesh(@test_mesh_2,opts);
print -dpng 'test_mesh_2.png';

display('Test mesh 3');
plot_lf_mesh(@test_mesh_3,opts);
print -dpng 'test_mesh_3.png';

display('Test mesh 4');
plot_lf_mesh(@test_mesh_4,opts);
print -dpng 'test_mesh_4.png';

display('Test mesh 5');
plot_lf_mesh(@test_mesh_5,opts);
print -dpng 'test_mesh_5.png';



