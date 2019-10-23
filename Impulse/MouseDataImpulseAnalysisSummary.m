%% Mouse Data Impulse Analysis Summary
% This function takes in the output of MouseDataImpulseAnalysis, a struct 
% called ImpulseAnalyzed that contains the traces that were analyzed as
% well as the linear fits and parameters.

% This function displays the summary table and summary graph for the user
% to see.

%Made by Andrianna Ayiotis
%Last updated 07/02/2018
function MouseDataImpulseAnalysisSummary(ImpulseAnalyzed,show_plots)
disp(ImpulseAnalyzed.summary)
t = ImpulseAnalyzed.t;
Ll = ImpulseAnalyzed.Ll;
Lr = ImpulseAnalyzed.Lr;
Rl = ImpulseAnalyzed.Rl;
Rr = ImpulseAnalyzed.Rr;
maxvel = ImpulseAnalyzed.info.maxvel;
mouse = ImpulseAnalyzed.info.mouse;
if(show_plots)
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    plot(t,Ll.chair_clean,'k')
    hold on
    plot(t,Rl.chair_clean,'k')
    plot(t,Ll.chairfit,'k','LineWidth',2)
    plot(t,Rl.chairfit,'k','LineWidth',2)
    plot(t,Ll.eye_clean,'r')
    plot(t,Ll.eyefit,'r','LineWidth',2)
    plot(t,Rl.eye_clean,'m')
    plot(t,Rl.eyefit,'m','LineWidth',2)
    title('Leftward Eye Movements')
    xlabel('Time (s)')
    ylabel('Angular Velocity (dps)')
    axis([0.1 0.3 -0.5*maxvel 1.3*maxvel])
    hold off
    subplot(1,2,2)
    plot(t,Lr.chair_clean,'k')
    hold on
    p1 = plot(t,Rr.chair_clean,'k');
    plot(t,Lr.chairfit,'k','LineWidth',2)
    p2 = plot(t,Rr.chairfit,'k','LineWidth',2);
    p3 = plot(t,Lr.eye_clean,'r');
    p4 = plot(t,Lr.eyefit,'r','LineWidth',2);
    p5 = plot(t,Rr.eye_clean,'m');
    p6 = plot(t,Rr.eyefit,'m','LineWidth',2);
    title('Rightward Eye Movements')
    xlabel('Time (s)')
    ylabel('Angular Velocity (dps)')
    legend([p1(1),p2,p3(1),p4,p5(1),p6],{'Inverted Chair','Chair Fit','Left Eye','Left Eye Fit','Right Eye','Right Eye Fit'})
    axis([0.1 0.3 -1.3*maxvel 0.5*maxvel])
    hold off
    suptitle([mouse,' Yaw Impulse Response'])
end
end