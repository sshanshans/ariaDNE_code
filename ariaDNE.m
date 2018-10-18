function [H] = ariaDNE(meshname, bandwidth, Options)

% This function computes the ariaDNE value of a mesh surface.
% ariaDNE is a robustly implemented algorithm for Dirichlet Normal
% Energy, which measures how much a surface deviates from a plane.

% Input:
%       meshname     - the mesh .ply file
%       bandwidth    - the episilon value in the paper, which indicates
%                      the size of local influence in the weight function
%       Options      - distInfo:  'Geodeisic' (default) or 'Euclidean'
%                    - distance:  precomputed distance matrix
%                    - cutThresh: cut-off threshold for very small values
%                                 default is set to be 0, i.e. no values
%                                 will be ignored.

% Output:
%       H.normals    - approximated normals for each vertex
%       H.curvature  - approxiamted curvature for each vertex
%       H.localDNE   - local energy values for each vertex
%       H.dne        - ARIADNE value for the surface

% Author:
%       Shan Shan
%       Department of Mathematics
%       sshan@math.duke.edu
%       June 14, 2018

% default options
H.Opts.distInfo = 'Geodeisic';
H.Opts.distance = [];
H.Opts.bandwidth = bandwidth;
H.Opts.cutThresh = 0;

if(nargin < 2)
    bandwidth = [];
    Options = struct();
end

if(isempty(bandwidth))
    bandwidth = 0.08;
end

% load user specified options
fn = fieldnames(Options);
for j = 1 : length(fn)
    name = fn{j};
    value = getfield(Options,name);
    
    if       strcmpi(name,'distance')      H.Opts.distance = value;
    elseif   strcmpi(name,'distInfo')      H.Opts.distInfo = value;
    elseif   strcmpi(name,'cutThresh')     H.Opts.cutThresh = value;
    else     fprintf('ARIADNE.m: invalid options "%s" ignored. \n', name);
    end
end

% mesh loading: convert .ply to .mat files
G = Mesh('ply', meshname);

% mesh cleaning: remove unreferenced vertices, zero-area faces, and
% isolated vertices.
G.remove_unref_verts;
G.remove_zero_area_faces;
G.DeleteIsolatedVertex;

% mesh preparation: centeralize, normalize the mesh to have surface area 1,
% compute an initial estimate of normals 
Centralize(G,'ScaleArea');
[~, face_area] = ComputeSurfaceArea(G);
vert_area = (face_area'*G.F2V)/3;
[vnorm, ~] = ComputeNormal(G);
if size(vnorm,1) < size(vnorm,2)
    vnorm = vnorm';
end

points = G.V';
numPoints = G.nV; 
normals = zeros(numPoints,3);
curvature = zeros(numPoints,1);
curvature_nn = zeros(numPoints,1);

% compute or load pairwise distnace
if isempty(H.Opts.distance) 
    if strcmpi(H.Opts.distInfo, 'Geodeisic') 
        d_dist = graphallshortestpaths(Triangulation2AdjacencyWeighted(G));
    elseif strcmpi(H.Opts.distInfo, 'Euclidean')
        d_dist = squareform(pdist(points));
    else
        fprintf('ARIADNE.m: invalid distInfo options "%s" ignored. \n', H.Opts.distInfo);
    end
else 
    d_dist = H.Opts.distance;
end

% define the weight matrix
K = exp(-d_dist.^2/(bandwidth^2));

% for each vertex in the mesh, estimate its curvature via PCA
for jj = 1:numPoints
    neighbour = find(K(jj,:) > H.Opts.cutThresh);
    numNeighbours = length(neighbour);
    if numNeighbours <= 3
     fprintf('ARIADNE.m: Too few neighbor on vertex %d. \n', jj);
    end
    p = repmat(points(jj,1:3),numNeighbours,1) - points(neighbour,1:3);
    w = K(jj, neighbour);
    
    % build covariance matrix for PCA
    C = zeros(1,6);
    C(1) = sum(p(:,1).*(w').*p(:,1),1);
    C(2) = sum(p(:,1).*(w').*p(:,2),1);
    C(3) = sum(p(:,1).*(w').*p(:,3),1);
    C(4) = sum(p(:,2).*(w').*p(:,2),1);
    C(5) = sum(p(:,2).*(w').*p(:,3),1);
    C(6) = sum(p(:,3).*(w').*p(:,3),1);
    C = C ./ sum(w);
    
    Cmat = [C(1) C(2) C(3);...
        C(2) C(4) C(5);...
        C(3) C(5) C(6)];  
    
    % compute its eigenvalues and eigenvectors
    [v,d] = eig(Cmat);
    d = diag(d);
    
    % find the eigenvector that is closest to the vertex normal 
    v_aug = [v -v];
    diff = v_aug - repmat(vnorm(jj,:)',1,6);
    q = sum(diff.^2,1);
    [~,k] = min(q);
    
    % use that eigenvector to give an updated estimate to the vertex normal
    normals(jj,:) = v_aug(:,k)';
    k = mod(k-1,3)+1;
    
    % use the eigenvalue of that egienvector to estimate the curvature
    lambda = d(k);
    curvature(jj) = lambda/sum(d);
    curvature_nn(jj) = nnz(w>max(w)*1e-4);
end

% save the outputs
H.curvature = curvature;
H.normals = normals;
H.dne = sum(curvature.*vert_area');
H.localDNE = curvature.*vert_area';
end

