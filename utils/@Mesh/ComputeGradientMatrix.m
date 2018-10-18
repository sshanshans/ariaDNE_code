function [gradx,grady,gradz,Fn] = ComputeGradientMatrix(G)
V=G.V; %vertices
F=G.F; %faces
Nv=size(V,2); % number of vertices
Nf=size(F,2); %number of faces
Fn=G.ComputeFaceNormals; %face normals
if size(F,1)~=3
    error('Not a triangular mesh!');
end

grad=zeros(3,3*Nf);
for i=1:Nf
    e=V(:,F([2,3,1],i))-V(:,F([1,2,3],i)); % edges of a triangle [e1,e2,e3]
    Fni=Fn(:,i);
    A=-det([e,Fni; 1 1 1 1])/6; % area of triangle
    if A<1e-7
        warning('Tiny triangle detected. Ignoring triangle')
        A=1e-7;
    end
    grad(:,3*i-2:3*i)=[0,-Fni(3), Fni(2);Fni(3),0,-Fni(1);-Fni(2),Fni(1),0]*e/(2*A);
end

gradx=sparse(kron(1:Nf,[1,1,1]),F(:),grad(1,:),Nf,Nv);
grady=sparse(kron(1:Nf,[1,1,1]),F(:),grad(2,:),Nf,Nv);
gradz=sparse(kron(1:Nf,[1,1,1]),F(:),grad(3,:),Nf,Nv);