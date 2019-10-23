%% Make Path
%This small script adds all folders and subfolders at current level into
%the path.

function MakePath
    addpath(genpath(pwd))
end