%% Clean Data
%This script calls the MouseDataSineCleanData file to batch clean all
%of the data. Cleaning refers to the process of desaccading and removing
%erroneous cycles with blinks 

%Written by Andrianna Ayiotis
%Last updated 07/31/2018
%Last run on 07/31/2018 

function CleanData
%% Find names of all files that need to be cleaned
%Assuming you're in the folder with the script (Sinusoids)
cd Segmented
files=dir(fullfile(cd,'*.mat'));
cd ../
%Make folder if it doesn't exist
if(exist('Cleaned','dir')~=7) 
    mkdir Cleaned
    MakePath;
end
cd Cleaned
fname = {files(:).name}';
%% Clean all the files at once
for i = 90:length(fname)
    load(fname{i},'SegDat')
    CleanDat = MouseDataSineCleanData(SegDat);
    save(CleanDat.info.fname,'CleanDat')
end
cd ../
end