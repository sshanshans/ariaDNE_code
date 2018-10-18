function [V, F] = remove_unref_verts(G)
% REMOVE_UNREF_VERTS - Remove unreferenced vertices from mesh
%
% Method written by Julie Winchester (julie.winchester@duke.edu)

	nVert = size(G.V, 2);
	unrefVertIdx = [];
	for i = 1:nVert;
		if length(find(ismember(G.F,i))) == 0
			unrefVertIdx = [unrefVertIdx i];
		end
	end

	disp(unrefVertIdx);
	G.V(:, unrefVertIdx) = [];
	flipUnrefVertIdx = fliplr(unrefVertIdx);

	for j = 1:length(flipUnrefVertIdx)
		G.F(ismember(G.F,[flipUnrefVertIdx(j) : nVert])) = ...
			G.F(ismember(G.F,[flipUnrefVertIdx(j) : nVert])) - 1;	
	end

    G.F2V = G.ComputeF2V;
    G.V2E = G.ComputeV2E;
    G.nV = size(G.V,2);
    G.nF = size(G.F,2);
    G.nE = size(G.V2E,2);

end

