function [path,vlist,plist] = ComputeGeodesic(G, D, Vind, options)

% compute_geodesic_mesh - extract a discrete geodesic on a mesh
%
%   [path,vlist,plist] = compute_geodesic_mesh(D, vertex, face, x, options);
%
%   D is the set of geodesic distances.
%
%   path is a 3D curve that is the shortest path starting at x.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

%%% gradient descent on edges
% precompute the adjacency datasets
E2F = G.ComputeE2F;
A = G.ComputeAdjacency;

% initialize the paths
[w,f] = vertex_stepping(G.F, Vind, E2F, A, D);

vlist = [Vind;w]; 
plist = [1];
Dprev = D(Vind);

while true;
    % current triangle
    i = vlist(1,end);
    j = vlist(2,end);
    k = get_vertex_face(G.F(:,f),i,j);
    a = D(i); b = D(j); c = D(k);
    % adjacent faces
    f1 = get_face_face(E2F, f, i,k);
    f2 = get_face_face(E2F, f, j,k);
    % compute gradient in local coordinates
    Vind = plist(end); y = 1-Vind;
    gr = [a-c;b-c];
    % twist the gradient
    u = G.V(:,i) - G.V(:,k);
    v = G.V(:,j) - G.V(:,k);
    A = [u v]; A = (A'*A)^(-1);
    gr = A*gr;
    nx = gr(1); ny = gr(2);
    % compute intersection point
    etas = -y/ny;
    etat = -Vind/nx;
    s = Vind + etas*nx;
    t = y + etat*ny;
    if etas<0 && s>=0 && s<=1 && f1>0
        %%% CASE 1 %%%
        plist(end+1) = s;
        vlist(:,end+1) = [i k];
        % next face
        f = f1;
    elseif etat<0 && t>=0 && t<=1 && f2>0
        %%% CASE 2 %%%
        plist(end+1) = t;
        vlist(:,end+1) = [j k];
        % next face
        f = f2;
    else
        %%% CASE 3 %%%
        if a<=b
            z = i;            
        else
            z = j;
        end
        [w,f] = vertex_stepping( face, z, e2f, v2v, D);
        vlist(:,end+1) = [z w];  
        plist(end+1) = 1;   
    end
    Dnew = D(vlist(1,end))*plist(end) + D(vlist(2,end))*(1-plist(end));
    if Dnew==0 || (Dprev==Dnew && length(plist)>1)
        break;
    end
    Dprev=Dnew;
end

v1 = G.V(:,vlist(1,:));
v2 = G.V(:,vlist(2,:));

path = v1.*repmat(plist, [3 1]) + v2.*repmat(1-plist, [3 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,f] = vertex_stepping(F, Vind, E2F, A, D)

% adjacent vertex with minimum distance
v2v=find(A(Vind,:));
[~,I] = min( D(v2v{v}) ); w = v2v{v}(I);
f1 = E2F(v,w);
f2 = E2F(w,v);
if f1<0
    f = f2; return;
end
if f2<0
    f = f1; return;
end
z1 = get_vertex_face(F(:,f1),v,w);
z2 = get_vertex_face(F(:,f2),v,w);
if D(z1)<D(z2);
    f = f1; 
else
    f = f2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = get_vertex_face(f,v1,v2)

if nargin==2
    v2 = v1(2); v1 = v1(1);
end
k = setdiff(f, [v1 v2]);
if length(k)~=1
    error('Error in get_vertex_face');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = get_face_face(e2f, f, i,j)

f1 = e2f(i,j); f2 = e2f(j,i);
if f==f1
    f = f2;
else
    f = f1;
end



