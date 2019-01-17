function mystlPlot(v, f, t)
%STLPLOT is an easy way to plot an STL object
%V is the Nx3 array of vertices
%F is the Mx3 array of faces
%NAME is the name of the object, that will be displayed as a title

figure  % <--- commentami 
object.vertices = v;
object.faces = f;
patch(object,'FaceColor',       [0.8 0.8 1], ...
         'EdgeColor',       'k',        ... 
         'FaceLighting',    'none',     ...
         'AmbientStrength', 0.01,...
         'FaceAlpha', t);hold on

% Add a camera light, and tone down the specular highlighting
% camlight('headlight');
% material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);
grid on;

