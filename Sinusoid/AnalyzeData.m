%% Analyze Data
%This script calls the MouseDataSineAnalysis file to batch analyze all
%of the data. 

%Written by Andrianna Ayiotis
%Last updated 07/02/2018
%Last run on 08/02/2018 
function AnalyzeData
%% Find names of all files that need to be analyzed
%Assuming you're in the folder with the script (Sinusoids)
cd Cleaned
files=dir(fullfile(cd,'*.mat'));
cd ../
%Make folder if it doesn't exist
if(exist('Analyzed','dir')~=7)
    mkdir Analyzed
    MakePath;
end
cd Analyzed
fname = {files(:).name}';
%% Analyze all the files at once
for i = 1:length(fname)  
    load(fname{i},'CleanDat')
    SineAnalyzed = MouseDataSineAnalysis(CleanDat);
    save(SineAnalyzed.info.fname,'SineAnalyzed')
end
cd ../
end