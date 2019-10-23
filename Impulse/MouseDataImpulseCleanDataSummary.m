%% Mouse Data Impulse Clean Data Summary
% This function takes in the output of MouseDataImpulseCleanData, a struct 
% called CleanImpulseData that contains the clean traces that will be
% analyzed.

%Made by Andrianna Ayiotis
%Last updated 07/16/2018

function MouseDataImpulseCleanDataSummary(CleanImpulseData)
t = CleanImpulseData.t;
Ll = CleanImpulseData.Ll;
Lr = CleanImpulseData.Lr;
Rl = CleanImpulseData.Rl;
Rr = CleanImpulseData.Rr;
maxvel = CleanImpulseData.info.maxvel;
mouse = CleanImpulseData.info.mouse;
%% Make summary figures
subplot(2,1,1)
hold on
p1 = plot(t,Ll.chair_clean,'k');
plot(t,Rl.chair_clean,'k')
p2 = plot(t,Ll.eye_clean,'r');
p3 = plot(t,Rl.eye_clean,'m');
legend([p1(1),p2(1),p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
title('Leftward Eye Movements')
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
axis([t(1) t(end) -0.5*maxvel 1.3*maxvel])
hold off
subplot(2,1,2)
p1 = plot(t,Lr.chair_clean,'k');
hold on
plot(t,Rr.chair_clean,'k');
p2 = plot(t,Lr.eye_clean,'r');
p3 = plot(t,Rr.eye_clean,'m');
legend([p1(1),p2(1),p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
title('Rightward Eye Movements')
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
axis([t(1) t(end) -1.3*maxvel 0.5*maxvel])
suptitle([mouse,' Impulse Response'])
hold off
end