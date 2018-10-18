function pathSetup(BaseDirectory)

% pathSetup adds the various directories used by the ARIADNE code
% to the current path.  pathSetup assumes the current directory is
% the base directory unless another directory is specified

fprintf('pathSetup.m: setting ariaDNE paths ... \n');

if nargin==0
  Prefix  = [pwd filesep];
 else
   Prefix  = [BaseDirectory filesep];
end;

% choose your utility code by changing this line
appendpath(([Prefix]));
appendpath(([Prefix 'utils']));

fprintf('pathSetup.m: disabling case sensitivity warning ... \n');
warning('off','MATLAB:dispatcher:InexactMatch');

function appendpath(string)

fprintf('\t%s\\ \n', string);
addpath(genpath(string));

return;
