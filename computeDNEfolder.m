function result = computeDNEfolder(foldername, bandwidth, Options)
    [fileNameList, suffix] = getFileNames(foldername);
    result = zeros(1, length(fileNameList));
    for i = 1:length(fileNameList)
        meshname = [foldername '/' fileNameList{i} suffix];
        H = ARIADNE(meshname, bandwidth, Options);
        result(i) = H.dne;
    end
end