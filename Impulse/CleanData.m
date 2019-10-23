%% Clean Data
%This script calls the MouseDataImpulseCleanData file to batch clean all
%of the data. 

%Written by Andrianna Ayiotis
%Last updated 07/02/2018
%Last run on 07/16/2018 
function CleanData
%% Find names of all files that need to be segmented
%Assuming you're in the folder with the script (Velocity Ramp)
cd EMAOutput
files=dir(fullfile(cd,'*.mat'));
cd ../
if(exist('Cleaned','dir')~=7) %If a folder called Cleaned already exists
    mkdir Cleaned  
    MakePath;
end
cd Cleaned
fname = {files(:).name}';
%% Clean them all
for i = 1:length(fname)
    load(fname{i},'ImpulseData')
    mouse = fname{i}(1:6);
    CleanImpulseData = MouseDataImpulseCleanData(mouse,fname{i},ImpulseData);
    save(CleanImpulseData.info.fname,'CleanImpulseData')
end
cd ../
end