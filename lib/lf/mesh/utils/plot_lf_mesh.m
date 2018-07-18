function plot_lf_mesh(mesh_data)
% Plots a mesh written by the write_matlab() mesh utility function
% of LehrFEM++

% Get geometry information for (affine) 2D mesh 
[x,y,TRI,QUAD,EDS] = mesh_data();

% Bounding box
bb = [min(x) max(x) min(y) max(y)];
dims = [bb(2)-bb(1) ; bb(4)-bb(3)];
bb = [bb(1)-0.1*dims(1) bb(2)+0.1*dims(1) ...
      bb(3)-0.1*dims(2) bb(4)+0.1*dims(2)];

figure('name','LehreFEM++ mesh'); axis(bb); hold on;

% plot vertices
plot(x,y,'r*');

% plot edges
num_edges = size(EDS,1);
for k=1:num_edges
    ed = [x(EDS(k,1)) x(EDS(k,2)) ; y(EDS(k,1)) y(EDS(k,2))];
    plot( ed(1,:), ed(2,:),'b-');
    baryc = sum(ed,2)/2;
    text(baryc(1),baryc(2),sprintf('%2i',k),'color','k');
end

% plot triangles with numbers
num_tri = size(TRI,1);
for  k=1:num_tri
    trc = [x(TRI(k,1)),x(TRI(k,2)),x(TRI(k,3)); ...
           y(TRI(k,1)),y(TRI(k,2)),y(TRI(k,3))];
    baryc = sum(trc,2)/3;
    bm = repmat(baryc,1,3);
    shrunk_tri = 0.9*(trc-bm)+bm;
    shrunk_tri = [shrunk_tri,shrunk_tri(:,1)];
    plot(shrunk_tri(1,:),shrunk_tri(2,:),'g-');
    text(baryc(1),baryc(2),sprintf('%2i',TRI(k,4)),'color','g');
    str = {'0','1','2'};
    text(shrunk_tri(1,1:3),shrunk_tri(2,1:3),str,'color','k','fontsize',6,'horizontalalignment','center');
end

% Plot quads
num_quad = size(QUAD,1);
for  k=1:num_quad
    trc = [x(QUAD(k,1)),x(QUAD(k,2)),x(QUAD(k,3)),x(QUAD(k,4)); ...
           y(QUAD(k,1)),y(QUAD(k,2)),y(QUAD(k,3)),y(QUAD(k,4)) ];
    baryc = sum(trc,2)/4;
    bm = repmat(baryc,1,4);
    shrunk_quad = 0.9*(trc-bm)+bm;
    shrunk_quad = [shrunk_quad,shrunk_quad(:,1)];
    plot(shrunk_quad(1,:),shrunk_quad(2,:),'m-');
    text(baryc(1),baryc(2),sprintf('%2i',QUAD(k,5)),'color','m');
end
