function touch(dir_string)
%TOUCH Summary of this function goes here
%   Detailed explanation goes here

if ~exist(dir_string, 'dir')
    mkdir(dir_string);
end

end

