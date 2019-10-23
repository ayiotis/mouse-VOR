%% Mouse Table Maker
%This script compiles and saves a table with the combination data for each 
%mouse so that it can be analyzed. It also saves a summary figure of all of 
%the mice into the Figures folder.

%Written by Andrianna Ayiotis
%Last updated 07/16/2018
%Last run on 07/16/2018 
%% Make the table
mice = {'WUm252','WUm255','WUm276','WUm277','WUm279','WUm282','WUm283','WUm284','WUm295','WUm296','WUm297'}';
type = {'double het','cko','cko','double het','double het','cko','cko','single het','double het','double het','cko'}';
mouse_nums = zeros(length(mice),7);
cd Figures
for i = 1:length(mice)
    fname = [mice{i},'ImpulseDataAnalyzed.mat'];
    load(fname,'ImpulseAnalyzed');
    MouseDataImpulseAnalysisSummary(ImpulseAnalyzed,true);
    savefig([mice{i},'.fig'])
    saveas(gcf,[mice{i},'.jpg'])
    close;
    mouse_nums(i,:) = ImpulseAnalyzed.summary{5,:};
end
labels = [{'mouseID','type'},ImpulseAnalyzed.summary.Properties.VariableNames];
allmousetab = [cell2table([mice,type]),array2table(mouse_nums)];
allmousetab.Properties.VariableNames = labels;
%save('MouseSummary.mat','allmousetab')
%% Make Summary Figures of the Types of Mice
het_tab = allmousetab(contains(allmousetab.type,'het'),:);
cko_tab = allmousetab(contains(allmousetab.type,'cko'),:);
het_mice = unique(het_tab.mouseID);
cko_mice = unique(cko_tab.mouseID);

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(het_mice)
    fname = [het_mice{i},'ImpulseDataAnalyzed.mat'];
    load(fname,'ImpulseAnalyzed');
    t = ImpulseAnalyzed.t;
    Ll = ImpulseAnalyzed.Ll;
    Lr = ImpulseAnalyzed.Lr;
    Rl = ImpulseAnalyzed.Rl;
    Rr = ImpulseAnalyzed.Rr;
    maxvel = ImpulseAnalyzed.info.maxvel;
    mouse = ImpulseAnalyzed.info.mouse;
    subplot(2,length(het_mice),i)
    hold on
    plot(t,Ll.chair_clean,'k')
    plot(t,Rl.chair_clean,'k')
    plot(t,Ll.chairfit,'k','LineWidth',2)
    plot(t,Rl.chairfit,'k','LineWidth',2)
    plot(t,Ll.eye_clean,'r')
    plot(t,Ll.eyefit,'r','LineWidth',2)
    plot(t,Rl.eye_clean,'m')
    plot(t,Rl.eyefit,'m','LineWidth',2)
    title({mouse, 'Leftward Eye Movements'})
    if(i==1)
        ylabel('Angular Velocity (dps)')
    end
    axis([0.1 0.3 -0.5*maxvel 1.3*maxvel])
    hold off
    subplot(2,length(het_mice),i+length(het_mice))
    hold on
    plot(t,Lr.chair_clean,'k');
    p1 = plot(t,Rr.chair_clean,'k');
    plot(t,Lr.chairfit,'k','LineWidth',2);
    p2 = plot(t,Rr.chairfit,'k','LineWidth',2);
    p3 = plot(t,Lr.eye_clean,'r');
    p4 = plot(t,Lr.eyefit,'r','LineWidth',2);
    p5 = plot(t,Rr.eye_clean,'m');
    p6 = plot(t,Rr.eyefit,'m','LineWidth',2);
    title('Rightward Eye Movements')
    xlabel('Time (s)')
    if(i==1)
        ylabel('Angular Velocity (dps)')
    end
    axis([0.1 0.3 -1.3*maxvel 0.5*maxvel])
    hold off
end
legend([p1(1),p2,p3(1),p4,p5(1),p6],{'Inverted Chair','Chair Fit','Left Eye','Left Eye Fit','Right Eye','Right Eye Fit'})
suptitle('Het Mice Yaw Impulse VOR Response')
savefig('AllHetMiceYawImpulse.fig')
saveas(gcf,'AllHetMiceYawImpulse.jpg')

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(cko_mice)
    fname = [cko_mice{i},'ImpulseDataAnalyzed.mat'];
    load(fname,'ImpulseAnalyzed');
    t = ImpulseAnalyzed.t;
    Ll = ImpulseAnalyzed.Ll;
    Lr = ImpulseAnalyzed.Lr;
    Rl = ImpulseAnalyzed.Rl;
    Rr = ImpulseAnalyzed.Rr;
    maxvel = ImpulseAnalyzed.info.maxvel;
    mid = floor(length(t)/4*3);
    mouse = ImpulseAnalyzed.info.mouse;
    subplot(2,length(cko_mice),i)
    hold on
    plot(t,Ll.chair_clean,'k')
    plot(t,Rl.chair_clean,'k')
    plot(t,Ll.chairfit,'k','LineWidth',2)
    plot(t,Rl.chairfit,'k','LineWidth',2)
    plot(t,Ll.eye_clean,'r')
    plot(t,Ll.eyefit,'r','LineWidth',2)
    plot(t,Rl.eye_clean,'m')
    plot(t,Rl.eyefit,'m','LineWidth',2)
    title({mouse, 'Leftward Eye Movements'})
    if(i==1)
        ylabel('Angular Velocity (dps)')
    end
    axis([0.1 0.3 -0.5*maxvel 1.3*maxvel])
    hold off
    subplot(2,length(cko_mice),i+length(cko_mice))
    hold on
    plot(t,Lr.chair_clean,'k');
    p1 = plot(t,Rr.chair_clean,'k');
    plot(t,Lr.chairfit,'k','LineWidth',2);
    p2 = plot(t,Rr.chairfit,'k','LineWidth',2);
    p3 = plot(t,Lr.eye_clean,'r');
    p4 = plot(t,Lr.eyefit,'r','LineWidth',2);
    p5 = plot(t,Rr.eye_clean,'m');
    p6 = plot(t,Rr.eyefit,'m','LineWidth',2);
    title('Rightward Eye Movements')
    xlabel('Time (s)')
    if(i==1)
        ylabel('Angular Velocity (dps)')
    end
    axis([0.1 0.3 -1.3*maxvel 0.5*maxvel])
    hold off
end
legend([p1(1),p2,p3(1),p4,p5(1),p6],{'Inverted Chair','Chair Fit','Left Eye','Left Eye Fit','Right Eye','Right Eye Fit'})
suptitle('Cko Mice Yaw Impulse VOR Response')
savefig('AllCkoMiceYawImpulse.fig')
saveas(gcf,'AllCkoMiceYawImpulse.jpg')
cd ../