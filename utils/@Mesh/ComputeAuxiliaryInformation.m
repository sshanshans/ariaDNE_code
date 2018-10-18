function [newG,Aux] = ComputeAuxiliaryInformation(G,options)
%Wrapper function to compute relevant quantitative auxilliary quantities. Needed in
%order to do in parallel fashion. This will set all auxiliary data except
%for the name.
 
%% Get aux field if already exists, should only be if observer landmarks were computed
if isfield(G,'Aux')
    Aux = G.Aux;
else
    Aux = {};
end
 
%centralize
[Aux.Area,Aux.Center] = G.Centralize('ScaleArea');
 
%uniformizae and compute conformal factor
[Aux.UniformizationV,~,Aux.VertArea] = G.ComputeMidEdgeUniformization();
 
Aux = G.ComputeCurvatureFeatures(Aux,options);
Aux.LB = G.ComputeCotanLaplacian;
 
minds = [Aux.GaussMaxInds;Aux.GaussMinInds;Aux.ConfMaxInds];
minds = unique(minds);
Aux.DensityPnts = G.GeodesicFarthestPointSampling(1000,minds);
G.Aux = Aux;
newG = G;
end