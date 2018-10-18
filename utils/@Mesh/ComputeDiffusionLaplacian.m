function [W,D,nL] = ComputeDiffusionLaplacian(G, NN, bandwidth)
%COMPUTEDIFFUSIONLAPLACIAN: compute the diffusion Laplacian on a mesh
%   connect each vertex to its NN nearest neighbors, then divide all
%   distances by bandwidth

nV = max(size(G.V));

if nargin < 2
    NN = nV-1;
    bandwidth = 'autoTune';
end
if nargin < 3
    bandwidth = 'autoTune';
end

atria = nn_prepare(G.V');
[Indx, Dist] = nn_search(G.V',atria,(1:nV)',NN+1,-1,0.0);
Indx(:,1) = [];
Dist(:,1) = [];

if ~isnumeric(bandwidth)
    threshold = min(10,NN);
    autoTune = zeros(size(Dist));
    for j=1:size(autoTune,1)
        autoTune(j,:) = Dist(j,threshold)*Dist(Indx(j,:),threshold);
    end
    autoTune = sqrt(autoTune);
    W = sparse(repmat((1:nV)',NN,1), Indx(:), exp(-Dist.^2)./autoTune);
else
    W = sparse(repmat((1:nV)',NN,1), Indx(:), exp(-Dist.^2/bandwidth));
end

W = min(W,W');
D = sparse(1:nV, 1:nV, sum(W,2), nV, nV, nV);

invSqrtD = sparse(1:nV, 1:nV, 1./sqrt(sum(W,2)), nV, nV, nV);
nL = invSqrtD*W*invSqrtD;

end

