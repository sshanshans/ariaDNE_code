function [BV,BE, BF] = FindBoundaries(G)
% find the boundry vertices of G
% out: BV: indices of the boundary vertices
%      BE: Vx2 matrix containing the indices to the two neighboring
%      vertices of each vertex
% if ~isempty(G.BE) && ~isempty(G.BV)
%     BE=G.BE;
%     BV=G.BV;
%     return;
% end
%if isempty(G.E)
    [~,E]=G.ComputeAdjacencyMatrix;
%else
%    E=G.E;
%end
Nv = size(G.V,2);
[I,J] = find(E);
BV = sparse(Nv,1);
BE = sparse(Nv,2);
BF = [];
for m=1:length(I)
    i = I(m); j = J(m);
    if xor(E(i,j),E(j,i))
        loc = find(BE(i,:) == 0);
        if ~isempty(loc)
            BE(i,loc(1)) = j;
        end
        loc = find(BE(j,:) == 0);
        if ~isempty(loc)
            BE(j,loc(1)) = i;
        end
        BV(i) = 1;
        BV(j) = 1;
        BF = [BF E(i,j) E(j,i)];
    end
end
temp = BF(BF>0);
BF = sort(unique(temp));
BV=find(BV);