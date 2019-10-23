%% Mouse Data Processing
%This script allows you to process the data from its form as a .mat table
%with eye velocities to characterized cycle averages.
%To start there should already be a folder titled "EMAOutput" with .mat 
%files that each have an eye velocity table named SineData with columns 
%Time, Lz, Rz, and Chair.
%% Segment Data
SegmentData;
%% Dessacade and Select Cycles
CleanData;
%% Extract Sine Parameters
AnalyzeData;