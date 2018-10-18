function [normal,normalf] = ComputeNormal(G,vertex,face)

% compute_normal - compute the normal of a triangulation
%
%   [normal,normalf] = compute_normal(vertex,face);
%
%   normal(i,:) is the normal at vertex i.
%   normalf(j,:) is the normal at face j.
%
%   Copyright (c) 2004 Gabriel Peyr�

if nargin<3
    vertex = G.V;
    face = G.F;
end

%[vertex,face] = G.CheckFaceVertex(vertex,face);

nface = size(face,2);
nvert = size(vertex,2);
% normal = zeros(3, nvert);

% unit normals to the faces
normalf = crossp( vertex(:,face(2,:))-vertex(:,face(1,:)), ...
                  vertex(:,face(3,:))-vertex(:,face(1,:)) );
d = sqrt( sum(normalf.^2,1) ); 
d(d<eps)=1;
normalf = normalf ./ repmat( d, 3,1 );


[~, face_area] = ComputeSurfaceArea(G);

% unit normal to the vertex
normal = zeros(3,nvert);
for i=1:nface
    f = face(:,i);
    for j=1:3
%         a = mod(j,3) +1 ;
%         b = mod(j+1,3) +1 ;
%         vec_a = vertex(:,f(a)) - vertex(:,f(j));
%         vec_b = vertex(:,f(b)) - vertex(:,f(j));
%         theta = atan2(norm(cross(vec_a',vec_b')),dot(vec_a',vec_b')) + pi;
        normal(:,f(j)) = normal(:,f(j)) + normalf(:,i); %* face_area(i)*theta; % SS modification: add weights of surface area
    end
end
% normalize
d = sqrt( sum(normal.^2,1) ); d(d<eps)=1;
normal = normal ./ repmat( d, 3,1 );

% enforce that the normal are outward
v = vertex - repmat(mean(vertex,1), 3,1);
s = sum( v.*normal, 2 );
if sum(s>0)<sum(s<0)
    % flip
    normal = -normal;
    normalf = -normalf;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(1,:) = x(2,:).*y(3,:) - x(3,:).*y(2,:);
z(2,:) = x(3,:).*y(1,:) - x(1,:).*y(3,:);
z(3,:) = x(1,:).*y(2,:) - x(2,:).*y(1,:);
