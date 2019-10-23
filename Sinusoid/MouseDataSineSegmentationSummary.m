%% Mouse Data Sine Segmentation Summary
%This function takes in the struct created by Mouse Data Sine Segmentation
%and displays the data so the user knows what to expect from that file.

%Written by Andrianna Ayiotis
%Last updated 06/13/2018

function MouseDataSineSegmentationSummary(SegDat)
%All the info the user could want
disp(SegDat.info)
%Plot the Data
round_freq = SegDat.info.round_freq;
mouse = SegDat.info.mouse;
maxvel = SegDat.info.maxvel;
plot_title = [num2str(mouse),' ',num2str(round_freq),' Hz Segmented Data'];
figure;
p1 = plot(SegDat.t,SegDat.Chair,'k');
hold on
p2 = plot(SegDat.t,SegDat.LEye,'r');
p3 = plot(SegDat.t,SegDat.REye,'m');
hold off
title(plot_title)
xlabel('Time (s)')
ylabel('Velocity (dps)')
axis([0 SegDat.t(end) -1.1*maxvel 1.1*maxvel]) 
legend([p1(1),p2(1),p3(1)],{'Chair','Left Eye','Right Eye'})
end