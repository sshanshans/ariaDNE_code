% set-up path:
clear; clc;
%%
pathSetup();
% pathSetup(BaseDirectory) %or provide a specified base directory

% load mesh
meshname = 'data2.ply';
G = Mesh('ply', meshname);G.remove_unref_verts;
G.remove_zero_area_faces;
G.DeleteIsolatedVertex;

% compute ARIADNE
Options.distInfo = 'Geodeisic';
Options.cutThresh = 0;
bandwidth = 0.04;
H = ARIADNE(meshname, bandwidth, Options);
fprintf('ARIADNE for data.ply is %f. \n', H.dne);

figure;
patch('Faces', G.F', 'Vertices', G.V', 'FaceVertexCData', H.localDNE, 'FaceColor','interp');
axis off;
axis equal;
cameratoolbar;
caxis([min(H.localDNE) 0.2*max(H.localDNE)]);


