%% Mouse Data Sine Analysis Summary
%This function takes in the struct created by Mouse Data Sine Analysis
%and displays the graphs and parameters made by the analysis.

%Written by Andrianna Ayiotis
%Last updated 07/31/2018
function MouseDataSineAnalysisSummary(SineAnalyzed,show_plots)
disp(SineAnalyzed.info)
disp(SineAnalyzed.sum_tab)
L = SineAnalyzed.L;
R = SineAnalyzed.R;
Chair = SineAnalyzed.Chair;
t = SineAnalyzed.t; 
mouse = SineAnalyzed.info.mouse;
freq = SineAnalyzed.info.round_freq;
maxvel = SineAnalyzed.info.CleanDat.maxvel;
if(show_plots)
    hold on    
    p1 = plot(t,Chair.cycle_avg,'k');
    plot(t,L.All,'Color',[1,0.7,0.7]) %color is light red
    p2 = plot(t,L.Fit,'Color','r');
    p3 = plot(t,L.eyefit,'r','LineWidth',3);
    plot(t,R.All,'Color',[1,0.7,1]) %color is light pink
    p4 = plot(t,R.Fit,'Color','m'); 
    p5 = plot(t,R.eyefit,'m','LineWidth',3);
    legend([p1,p2(1),p3,p4(1),p5],{'Chair','Left Traces','Left Fit','Right Traces','Right Fit'})
    suptitle([mouse,' ',num2str(freq),'Hz'])
    axis([t(1),t(end),-maxvel,maxvel])
    hold off
end
end