% This is an example computing ariaDNE for 'data.ply'

% set-up path:
clear; clc;
pathSetup();
% pathSetup(BaseDirectory) %or provide a specified base directory

% compute ariaDNE
Options.distInfo = 'Geodeisic';
Options.cutThresh = 0;
meshname = 'data.ply';
bandwidth = 0.08;
H = ariaDNE(meshname, bandwidth, Options);
fprintf('ariaDNE for data.ply is %f. \n', H.dne);



