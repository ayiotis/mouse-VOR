%% Make Path
% This small function adds all folders and subfolders on your current path 
% to the path programatically.
% Helpful if you are making a new folder programmatically and want to add
% it to the path.
function MakePath
    addpath(genpath(pwd))
end