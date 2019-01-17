function stlPlot(v, f, name, viewparam)
%STLPLOT is an easy way to plot an STL object
%V is the Nx3 array of vertices
%F is the Mx3 array of faces
%NAME is the name of the object, that will be displayed as a title

figure;
object.vertices = v;
object.faces = f;
patch(object,'FaceColor',       [0.8 0.0 0.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
az = viewparam.az;
el = viewparam.el;
view(az, el);
grid off;
axis off;
title(name);
