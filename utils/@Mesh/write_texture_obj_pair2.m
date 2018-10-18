function write_texture_obj_pair(GM, GN, cpDistStructMN, cpDistStructNM, filePath)
% WRITE_TEXTURE_OBJ_PAIR - Save 2 .obj files with surface map as texture

	VM = GM.V;
	VN = GN.V;

	GM.V = GM.V*4;
	GN.V = GN.V*4;
	pathM = fullfile(filePath, [GM.Aux.name '_map.obj']);
	disp(pathM);
	pathN = fullfile(filePath, [GN.Aux.name '_map.obj']);
	disp(pathN);
	if cpDistStructMN.cPdist<cpDistStructNM.cPdist
	    optionsM.Texture.Coordinates = cpDistStructMN.TextureCoords1/2+0.5;
	    optionsN.Texture.Coordinates = cpDistStructMN.TextureCoords2/2+0.5;
	else
	    optionsM.Texture.Coordinates = cpDistStructNM.TextureCoords2/2+0.5;
	    optionsN.Texture.Coordinates = cpDistStructNM.TextureCoords1/2+0.5;
	end
	GM.Write(pathM,'obj',optionsM); 
	GN.Write(pathN,'obj',optionsN);

	GM.V = VM;
	GN.V = VN;

end