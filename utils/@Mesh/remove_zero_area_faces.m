function F = remove_zero_area_faces(G)
% REMOVE_ZERO_AREA_FACES - Given a  @Mesh object, return a face array with no zero area faces
%
% Method written by Julie Winchester (julie.winchester@duke.edu)

	zeroFaceInds = [];
	for i = 1:size(G.F, 2);
		faceVerts = G.V(:, G.F(:, i));
		A = triangle_area(faceVerts(:, 1), faceVerts(:, 2), faceVerts(:, 3));
		if A == 0
			disp(['Zero area face found at index ' num2str(i)]);
			zeroFaceInds = [zeroFaceInds i];
		end
	end
	disp(G.F(:, zeroFaceInds));
	G.F(:, zeroFaceInds) = [];

	function A = triangle_area(a, b, c)
		A = 0.5 * sqrt(sum(abs(cross(b-a, c-a)).^2));
	end

    G.F2V = G.ComputeF2V;
    G.V2E = G.ComputeV2E;
    G.nV = size(G.V,2);
    G.nF = size(G.F,2);
    G.nE = size(G.V2E,2);

end
