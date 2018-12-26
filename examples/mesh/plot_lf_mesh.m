function plot_lf_mesh(mesh_data,opts)
% Plots a mesh written by the write_matlab() mesh utility function
% of LehrFEM++
% mesh_data must be a function handle to a MATLAB function produced by 
% the writeMatlab() C++ function of the LehrFEM++ library.
% 
% The following options are available:
% opts.numbers -> Add entity index numbers
% opts.parents -> Add index numbers of parents
%   In this case opts.parents must be a handle to a function
%   produced by the LehrFEM++ function writeMatlabLevel()
%   (available in the refinement module)    
%

% Get geometry information for (affine) 2D mesh 
[x,y,TRI,QUAD,EDS] = mesh_data();

if (nargin < 2), opts = []; end

% Option: read parent information
if (isfield(opts,'parents'))
    display('Plotting mesh with parent information');
    [PTPAR,EDPAR,CELLPAR] = opts.parents();
    if (size(PTPAR,1) ~= length(x))
        error('ParentInfo Mismatch: number of points');
    end
    if (size(EDS,1) ~= size(EDPAR,1))
        error('ParentInfo Mismatch: number of edges');
    end
    if ((size(TRI,1) + size(QUAD,1)) ~= size(CELLPAR,1))
        error('ParentInfo Mismatch: number of cells')
    end
    parplot = true;
else
    PTPAR = []; EDPAR = []; CELLPAR = [];
    parplot = false;
end


% Bounding box
bb = [min(x) max(x) min(y) max(y)];
dims = [bb(2)-bb(1) ; bb(4)-bb(3)];
bb = [bb(1)-0.1*dims(1) bb(2)+0.1*dims(1) ...
      bb(3)-0.1*dims(2) bb(4)+0.1*dims(2)];

figure('name','LehreFEM++ mesh'); axis(bb); hold on;

% plot vertices
if (isfield(opts,'vertexmark'))
    plot(x,y,'r*');
else
    num_vertices = length(x);
    for k=1:num_vertices
        text(x(k),y(k),sprintf('%i',k-1),'Color','red',...
             'horizontalalignment','center','fontsize',8);
    end
end 

% plot edges
num_edges = size(EDS,1);
for k=1:num_edges
    ed = [x(EDS(k,1)) x(EDS(k,2)) ; y(EDS(k,1)) y(EDS(k,2))];
    dir = ed(:,2) - ed(:,1);
    quiver(ed(1,1),ed(2,1),dir(1),dir(2),'b-');
    %    plot( ed(1,:), ed(2,:),'b-');
    baryc = sum(ed,2)/2;
    if (isfield(opts,'numbers'))
        text(baryc(1),baryc(2),sprintf('%2i',k-1),...
             'fontsize',8,'color','b','horizontalalignment','center');
    end
    if (parplot)
        text(baryc(1),baryc(2),sprintf('P=%2i',EDPAR(k,1)),...
             'fontsize',8,'color','k','horizontalalignment','center');
    end        
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
    if (isfield(opts,'numbers'))
        text(baryc(1),baryc(2),sprintf('%2i',TRI(k,4)),...
             'fontsize',8,'color','g','horizontalalignment','center');
    end
    if (parplot)
        % indicate number of parent
        text(baryc(1),baryc(2),sprintf('P=%2i',CELLPAR(TRI(k,4)+1,1)),...
             'fontsize',8,'color','g','horizontalalignment','center');
        % mark refinement edge
        refedge = CELLPAR(TRI(k,4)+1,3);
        refedgepos = (shrunk_tri(:,refedge+1)+shrunk_tri(:,refedge+2))/2;
        plot(refedgepos(1),refedgepos(2),'r^');
    end    
    if (isfield(opts,'numbers'))
        str = {'0','1','2'};
        text(shrunk_tri(1,1:3),shrunk_tri(2,1:3),str,...
             'fontsize',8,'color','k','fontsize',6,'horizontalalignment','center');
    end
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
    if (isfield(opts,'numbers'))
        text(baryc(1),baryc(2),sprintf('%2i',QUAD(k,5)),...
             'fontsize',8,'color','m','horizontalalignment','center');
    end
    if (parplot)
        text(baryc(1),baryc(2),sprintf('P=%2i',CELLPAR(QUAD(k,5)+1,1)),...
             'fontsize',8,'color','m','horizontalalignment','center');
    end
    if (isfield(opts,'numbers'))
        str = {'0','1','2','3'};
        text(shrunk_quad(1,1:4),shrunk_quad(2,1:4),str,...
             'fontsize',8,'color','k','fontsize',6,'horizontalalignment','center');
    end
end
