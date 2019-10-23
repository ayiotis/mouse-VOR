%% Analyze Data
%This script calls the MouseDataImpulseAnalysis file to batch analyze all
%of the data. 
%Written by Andrianna Ayiotis
%Last updated 07/10/2018
%Last run on 07/10/2018 
function AnalyzeData
%% Find names of all files that need to be analyzed
%Assuming you're in the folder with the script (Velocity Ramp)
cd Cleaned
files=dir(fullfile(cd,'*.mat'));
cd ../
if(exist('Analyzed','dir')~=7)
    mkdir Analyzed
    MakePath;
end
cd Analyzed
fname = {files(:).name}';
%% Analyze all the files at once
plot_check = false;
for i = 1:length(fname)  
    load(fname{i},'CleanImpulseData')
    ImpulseAnalyzed = MouseDataImpulseAnalysis(CleanImpulseData,plot_check);
    save(ImpulseAnalyzed.info.fname,'ImpulseAnalyzed')
end
cd ../
end