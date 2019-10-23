%% Mouse Data Sine Clean Data Summary
%This function takes in the struct created by Mouse Data Sine Clean Data
%and displays the data so the user knows what to expect from that file.
%These data are ready for analysis.

%Written by Andrianna Ayiotis
%Last updated 06/21/2018

function MouseDataSineCleanDataSummary(CleanDat)
%All the info the use could want
disp(CleanDat.info)
%Plot the Data
round_freq = CleanDat.info.round_freq;
mouse = CleanDat.info.mouse;
t = CleanDat.t;
m_chair = CleanDat.Chair.cycle_avg;
Leye_seg = CleanDat.LEye.cycles;
Reye_seg = CleanDat.REye.cycles;
figure;
p1 = plot(t,m_chair,'k','LineWidth',2);
hold on
p2 = plot(t,Leye_seg,'Color','r'); 
p3 = plot(t,Reye_seg,'Color','m'); 
hold off
legend([p1,p2(1),p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
long_title = [mouse,' ',num2str(round_freq),' Hz Cycles'];
title(long_title)
xlabel('Time (s)')
ylabel('Velocity (dps)')
end