function [GPLmkIdx,ptuq] = GetGPLmk(G,numLmk)
%GETGPLMK Summary of this function goes here
%   Detailed explanation goes here

% if nargin < 3
%     lambda = 0.5;
% end
%normalize the mesh
[~,TriArea] = G.ComputeSurfaceArea();  %compute the area of each triangle
G.Aux.VertArea = G.F2V'*TriArea; %compute area of traingles adjacent to each vertex

% tingran wrote a different way to pick the curvature. 
% [Cgauss,Cmean] = G.ComputeCurvature();
% Lambda = G.Aux.VertArea.*(lambda*abs(Cgauss)/sum(abs(Cgauss))+(1-lambda)*abs(Cmean)/sum(abs(Cmean)));

%[vnorm, ~] = ComputeNormal(G);
%[curvature, ~, ~] = findPCANormalsImpr(G.V, vnorm ,0.03);
%curvature = exp(-G.V(3,:)');
%Lambda = curvature;
[~,curvature] = findPointNormals(G.V',10); %find the "curvature" using PCA method
Lambda = G.Aux.VertArea.*curvature/sum(curvature); %define curvature weight Lambda
%curvature = G.Aux.Conf;
%Lambda = log(G.Aux.VertArea.*curvature/10*sum(curvature));

% some auxillary stuffs that aren't so important
if size(G.E,1) > 2
    [I,J] = find(tril(G.E));
    G.E = ([I,J])';
end
EdgeIdxI = G.E(1,:);
EdgeIdxJ = G.E(2,:);
bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/5;
display(bandwidth);

% construct a sparse kernel for fast computation
BNN = min(500,G.nV); 
atria = nn_prepare(G.V');
[idx, dist] = nn_search(G.V',atria,(1:G.nV)',BNN+1,-1,0.0);
fullPhi = sparse(repmat(1:G.nV,1,BNN+1),idx,exp(-dist.^2/bandwidth),G.nV,G.nV);

% PDistMat = squareform(pdist(G.V'));
% fullPhi = exp(-PDistMat.^2/bandwidth);

% construct full kernel = Phi * Curvature Weight * Phi
disp('Constructing full kernel......');
tic;
fullMatProd = fullPhi * sparse(1:G.nV,1:G.nV,Lambda,G.nV,G.nV) * fullPhi;
disp(['full kernel constructed in ' num2str(toc) ' sec.']);

% kernel trace for initial landmark point 
KernelTrace = diag(fullMatProd);

%define the GPLmkIdx so that you can store the landmark Id later
GPLmkIdx = zeros(1,numLmk);
%same purpose as before 
invKn = zeros(numLmk);

%actual landmark picking process
cback = 0;
for j=1:numLmk
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d\n',j);
    
    %pick the first landmark: ptuq here gives the uncertatiny score for
    %each vertex
    if j == 1
        ptuq = KernelTrace;
    else
        %for the 2nd landmark
        if j == 2
            %define the inverse of the K(Tn, Tn);
            invKn(1:(j-1),1:(j-1)) = 1/fullMatProd(GPLmkIdx(1),GPLmkIdx(1));
            %calculate uncertainty score
            ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
                .*(invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:)),1)';
        else
            %for landmark >=3
            p = fullMatProd(GPLmkIdx(1:(j-2)),GPLmkIdx(j-1));
            mu = 1./(fullMatProd(GPLmkIdx(j-1),GPLmkIdx(j-1))-p'*invKn(1:(j-2),1:(j-2))*p);
            invKn(1:(j-2),1:(j-1)) = invKn(1:(j-2),1:(j-2))*[eye(j-2)+mu*(p*p')*invKn(1:(j-2),1:(j-2)),-mu*p];
            invKn(j-1,1:(j-1)) = [invKn(1:(j-2),j-1)',mu];
            productEntity = invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:);
            ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
                .*productEntity,1)';
        end
    end
    [~,maxUQIdx] = max(ptuq);%pick the landmark with the largest uncertainty 
    GPLmkIdx(j) = maxUQIdx; %store it into GPLmkIdx
end

%recompute the inv and ptuq
p = fullMatProd(GPLmkIdx(1:(end-1)),GPLmkIdx(end));
mu = 1./(fullMatProd(GPLmkIdx(end),GPLmkIdx(end))-p'*invKn(1:(end-1),1:(end-1))*p);
invKn(1:(end-1),:) = invKn(1:(end-1),1:(end-1))*[eye(numLmk-1)+mu*(p*p')*invKn(1:(end-1),1:(end-1)),-mu*p];
invKn(end,:) = [invKn(1:(end-1),end)',mu];
ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx)'...
    .*(invKn*fullMatProd(GPLmkIdx,:)),1)';

end

