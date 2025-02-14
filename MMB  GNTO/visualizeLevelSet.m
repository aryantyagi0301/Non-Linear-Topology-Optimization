% This function is from the paper "Zhang W, Li D, Yuan J, Song J, Guo X 
% (2016) A new three-dimensional topology optimization method based on 
% moving morphable components (MMCs). Comput Mech 59:647-665. 
% https://doi.org/10.1007/s00466-016-1365-0"
function h = visualizeLevelSet(g, data, displayType, level, titleStr)
% visualizeLevelSet: Display the level set at a particular time.
%
%   h = visualizeLevelSet(g, data, displayType, level, titleStr)
%
% Displays a variety of level set visualizations in dimensions 1 to 3.
%   The current figure and axis is used.
%
% A warning will be generated if the requested level set is missing.
%   For those display types that do not plot a level set
%   (such as a surf plot in 2D), the warning can be disabled by 
%   setting parameter level to be the empty vector [].
%
% Parameters:
%   g   	 Grid structure.
%   data         Array storing the implicit surface function.
%   displayType  String specifying the type of visualization (see below).
%   level        Which isosurface to display.  Defaults to 0.
%   titleStr     Optional string to place in the figure title.
%
%   h            Handle to the graphics object created.
%
% Display type options:
%
% Dimension 1:
%   'plot'       Plot the function value vs the state as a line.  If the 
%                  number of grid nodes is small, mark the node value as well.
%
% Dimension 2:
%   'contour'    Show an isocontour of the function as a solid curve in
%                  two dimensions.  Permits vector values for level parameter.
%   'surf'       Plot the function value vs the two dimensional state as
%                  a surface plot.
%
% Dimension 3:
%    'surface'   Show an isosurface of the function as a solid surface in
%                  three dimensions.
%    'slice'     On slices through the x,y,z midpoints of the grid show
%                  the function value through color coding.
%    'contourslice'  On slices through the x,y,z midpoints of the grid
%                  show an isocontour of the function as a solid curve.
%    'wireframe' Show an isosurface as a wireframe.  Essentially the same
%                  as surface, but with the faces turned off and the edges
%                  turned on.
%=========================================================
 [ mesh_xs, mesh_data ] = gridnd2mesh(g, data);
 h   = patch(isosurface(mesh_xs{:}, mesh_data, level));
 set( h, 'FaceColor', 'c', 'EdgeColor', 'none','facealpha',1);
 h2 = patch(isocaps(mesh_xs{:}, mesh_data, level));
 set(h2, 'FaceColor', 'c',	'EdgeColor','none');
 colormap([1 0 0])
 lighting phong;
 view(3)
 axis image;axis off;
 camlight right;
 pause(1e-6);
     
     