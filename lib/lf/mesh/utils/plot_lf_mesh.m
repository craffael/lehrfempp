function plot_lf_mesh(mesh_data)
% Plots a mesh written by the write_matlab() mesh utility function
% of LehrFEM++

% Get geometry information for (affine) 2D mesh 
[x,y,TRI,QUAD,EDS] = mesh_data();

figure('name','LehreFEM++ mesh');
axis([min(x) max(x) min(y) max(y)]);
hold on;

% plot vertices
plot(x,y,'r*');

% plot edges
num_edges = size(EDS,1);
for k=1:num_edges
    plot([x(EDS(k,1)) x(EDS(k,2))],[y(EDS(k,1)) y(EDS(k,2))],'b-');
end

