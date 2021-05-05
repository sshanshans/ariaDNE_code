function [Area,TriArea] = ComputeSurfaceArea(G)

V=G.V';
F=G.F';
L1 = sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2));
L2 = sqrt(sum((V(F(:,1),:)-V(F(:,3),:)).^2,2));
L3 = sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2));
S=(L1+L2+L3)/2;
temp = S.*(S-L1).*(S-L2).*(S-L3);
% uncomment the following line to get rid of really small triangles
%temp(temp<eps) = 0;
TriArea=sqrt(temp);
Area=sum(TriArea);

end
