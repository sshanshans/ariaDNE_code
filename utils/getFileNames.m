function [fileNameList, suffix] = getFileNames(pathName)
%GETFILES Summary of this function goes here
%   Detailed explanation goes here

tmp = dir(pathName);
fileNameList = arrayfun(@(x) x.name, tmp, 'UniformOutput', 0)';
fileNameList = fileNameList(3:end);
[fileNameList, suffix] = cellfun(@(x) strtok(x, '.'), fileNameList, 'UniformOutput', 0);
suffix = suffix{1};

end

