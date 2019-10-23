%% Mouse Data Processing
%This script allows you to process the data from its form as a .mat table
%with eye velocities to characterized cycle averages.
%To start there should already be a folder titled "EMAOutput" with .mat 
%files that each have an eye velocity table named ImpulseData with columns 
%Time, Lz, Rz, and Chair.
%% Dessacade and Select Cycles
CleanData;
%% Extract Sine Parameters
AnalyzeData;