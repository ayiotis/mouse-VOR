%% Mouse Stats
%This script analyzes the data in the table in MouseSummary.mat and saves a
%table with a summary of all the statistics run. It also saves relavant
%summary plots to the Figures folder.

%Written by Andrianna Ayiotis
%Last updated 07/20/2018
%Last run on 07/20/2018 
%% Load the table
load('MouseSummary.mat','allmousetab')
het_tab = allmousetab(contains(allmousetab.type,'het'),:);
cko_tab = allmousetab(contains(allmousetab.type,'cko'),:);
het_mice = unique(het_tab.mouseID);
cko_mice = unique(cko_tab.mouseID);

het_Ga = het_tab.meanGa;
cko_Ga = cko_tab.meanGa;
het_lat = het_tab.meanlat;
cko_lat = cko_tab.meanlat;
het_AMGa = het_tab.meanAmericoGa;
cko_AMGa = cko_tab.meanAmericoGa;

% Non-parametric test (Mann-Whitney U-test)
p1 = ranksum(het_Ga,cko_Ga);
p2 = ranksum(het_lat,cko_lat);
p3 = ranksum(het_AMGa,cko_AMGa);
% Parametric test (Two-sided t-test with unequal variance)
[~,p4] = ttest2(het_Ga,cko_Ga,'Vartype','unequal');
[~,p5] = ttest2(het_lat,cko_lat,'Vartype','unequal');
[~,p6] = ttest2(het_AMGa,cko_AMGa,'Vartype','unequal');
%% Without the outlier (WUm252)
het_tab2 = het_tab(~contains(het_tab.mouseID,'WUm252'),:);
het_Ga2 = het_tab2.meanGa;
het_lat2 = het_tab2.meanlat;
het_AMGa2 = het_tab2.meanAmericoGa;

% Non-parametric test (Mann-Whitney U-test)
[p7,~,stat7] = ranksum(het_Ga2,cko_Ga);
[p8,~,stat8] = ranksum(het_lat2,cko_lat);
[p9,~,stat9] = ranksum(het_AMGa2,cko_AMGa);

% Parametric test (Two-sided t-test with unequal variance)
[~,p10,~,stat10] = ttest2(het_Ga2,cko_Ga,'Vartype','unequal');
[~,p11,~,stat11] = ttest2(het_lat2,cko_lat,'Vartype','unequal');
[~,p12,~,stat12] = ttest2(het_AMGa2,cko_AMGa,'Vartype','unequal');
%% Compile all statistics
p = [p1,p2,p3;p4,p5,p6;p7,p8,p9;p10,p11,p12];
stat_sum = [cell2table({'Mann-Whitney U-test','All';'Two-Sided t-test','All';'Mann-Whitney U-test','NoOutliers';'Two-Sided t-test','NoOutliers'}),array2table(p)];
stat_sum.Properties.VariableNames = {'StatisticalTest','DataUsed','Ga','Latency','AMGa'};

disp('p-values:')
disp(stat_sum)
save('MouseStatisticsSummary.mat','stat_sum')
cd Figures
%% Make figures with the data points
mouse_groups = {'Het','Cko'};
spac = 0.01;
x = spac*((0:length(het_Ga)-1)-length(het_Ga)+1)+1-0.05;
x1 = spac*((0:length(cko_Ga)-1)-length(cko_Ga)+1)+1.2-0.05;

figure('units','normalized','outerposition',[0 0 0.5 1])
subplot(1,2,1)
scatter(x,het_Ga,[],'ko')
hold on
errorbar(x,het_Ga,het_tab.stdGa,'k.')
scatter(x1,cko_Ga,'b^')
errorbar(x1,cko_Ga,cko_tab.stdGa,'b.')
scatter(1,mean(het_Ga),'k','filled')
errorbar(1,mean(het_Ga),std(het_Ga),'k.')
scatter(1.2,mean(cko_Ga),'b^','filled')
errorbar(1.2,mean(cko_Ga),std(cko_Ga),'b.')
axis([0.875 1.225 0.8 1.6])
set(gca,'XTick',[1,1.2],'XTickLabel',mouse_groups)
xlabel('Mouse Type')
ylabel('Constant Acceleration Gain (Eye Acceleration/Head Acceleraion)')
title('Constant Acceleration Gain for Each Mouse')
hold off
subplot(1,2,2)
f1=scatter(x,het_lat,[],'k');
hold on
errorbar(x,het_lat,het_tab.stdlat,'k.')
f2=scatter(x1,cko_lat,'b^');
errorbar(x1,cko_lat,cko_tab.stdlat,'b.')
f3=scatter(1,mean(het_lat),'k','filled');
errorbar(1,mean(het_lat),std(het_lat),'k.')
f4=scatter(1.2,mean(cko_lat),'b^','filled');
errorbar(1.2,mean(cko_lat),std(cko_lat),'b.')
axis([0.875 1.225 0 55])
set(gca,'XTick',[1,1.2],'XTickLabel',mouse_groups)
xlabel('Mouse Type')
ylabel('Eye Response Latency (ms)')
title('Eye Response Latency for Each Mouse')
suptitle('Yaw Impulse Measurements') 
hold off
legend([f1,f3,f2,f4],{'Het Mice','Het Avg','Cko Mice','Cko Avg'})
savefig('ImpulseDataByGroup.fig')
saveas(gcf,'ImpulseDataByGroup.jpg')
%% Make Boxplots with the data points (without outlier)
x = spac*((0:length(het_Ga2)-1)-length(het_Ga2)+1)+1-0.05;
figure('units','normalized','outerposition',[0 0 0.5 1])
subplot(1,2,1)
scatter(x,het_Ga2,[],'ko')
hold on
errorbar(x,het_Ga2,het_tab2.stdGa,'k.')
scatter(x1,cko_Ga,'b^')
errorbar(x1,cko_Ga,cko_tab.stdGa,'b.')
scatter(1,mean(het_Ga2),'k','filled')
errorbar(1,mean(het_Ga2),std(het_Ga2),'k.')
scatter(1.2,mean(cko_Ga),'b^','filled')
errorbar(1.2,mean(cko_Ga),std(cko_Ga),'b.')
axis([0.875 1.225 0.8 1.6])
set(gca,'XTick',[1,1.2],'XTickLabel',mouse_groups)
xlabel('Mouse Type')
ylabel('Constant Acceleration Gain (Eye Acceleration/Head Acceleraion)')
title('Constant Acceleration Gain for Each Mouse')
hold off
subplot(1,2,2)
hold on
f1=scatter(x,het_lat2,[],'ko'); 
errorbar(x,het_lat2,het_tab2.stdlat,'k.');
f2=scatter(x1,cko_lat,'b^'); 
errorbar(x1,cko_lat,cko_tab.stdlat,'b.');
f3=scatter(1,mean(het_lat2),'k','filled');
errorbar(1,mean(het_lat2),std(het_lat2),'k.');
f4=scatter(1.2,mean(cko_lat),'b^','filled');
errorbar(1.2,mean(cko_lat),std(cko_lat),'b.');
axis([0.875 1.225 0 18])
legend([f1,f3,f2,f4],{'Het Mice','Het Avg','Cko Mice','Cko Avg'})
set(gca,'XTick',[1,1.2],'XTickLabel',mouse_groups)
xlabel('Mouse Type')
ylabel('Eye Response Latency (ms)')
title('Eye Response Latency for Each Mouse')
suptitle('Yaw Impulse Measurements without Outliers') 
hold off
savefig('ImpulseDataByGroupNoOutliers.fig')
saveas(gcf,'ImpulseDataByGroupNoOutliers.jpg')

rampfig.x = x;
rampfig.x1 = x1;
rampfig.het_Ga2 = het_Ga2;
rampfig.het_tab2 = het_tab2;
rampfig.cko_Ga = cko_Ga;
rampfig.cko_tab = cko_tab;
rampfig.mouse_groups = mouse_groups;
rampfig.het_lat2 = het_lat2;
rampfig.het_tab2 = het_tab2;
rampfig.cko_lat = cko_lat;
rampfig.cko_tab = cko_tab;
%% Make an AM Gain Comparrison 
x = spac*((0:length(het_AMGa)-1)-length(het_AMGa)+1)+1-0.05;
x1 = spac*((0:length(cko_AMGa)-1)-length(cko_AMGa)+1)+1.2-0.05;
x2 = spac*((0:length(het_AMGa2)-1)-length(het_AMGa2)+1)+1-0.05;

figure('units','normalized','outerposition',[0 0 0.5 1])
subplot(1,2,1)
scatter(x,het_AMGa,[],'k')
hold on
errorbar(x,het_AMGa,het_tab.stdAmericoGa,'k.')
scatter(x1,cko_AMGa,'b^')
errorbar(x1,cko_AMGa,cko_tab.stdAmericoGa,'b.')
scatter(1,mean(het_AMGa),'k','filled')
errorbar(1,mean(het_AMGa),std(het_AMGa),'k.')
scatter(1.2,mean(cko_AMGa),'b^','filled')
errorbar(1.2,mean(cko_AMGa),std(cko_AMGa),'b.')
axis([0.875 1.225 0.3 1.5])
set(gca,'XTick',[1,1.2],'XTickLabel',mouse_groups)
xlabel('Mouse Type')
ylabel('Constant Acceleration Gain (Eye Acceleration/Head Acceleraion)')
title('All Mice')
hold off
subplot(1,2,2)
f1=scatter(x2,het_AMGa2,[],'ko');
hold on
errorbar(x2,het_AMGa2,het_tab2.stdAmericoGa,'k.')
f2=scatter(x1,cko_AMGa,'b^');
errorbar(x1,cko_AMGa,cko_tab.stdAmericoGa,'b.')
f3=scatter(1,mean(het_AMGa2),'k','filled');
errorbar(1,mean(het_AMGa2),std(het_AMGa2),'k.')
f4=scatter(1.2,mean(cko_AMGa),'b^','filled');
errorbar(1.2,mean(cko_AMGa),std(cko_AMGa),'b.')
axis([0.875 1.225 0.3 1.5])
legend([f1,f3,f2,f4],{'Het Mice','Het Avg','Cko Mice','Cko Avg'})
set(gca,'XTick',[1,1.2],'XTickLabel',mouse_groups)
xlabel('Mouse Type')
ylabel('Constant Acceleration Gain (Eye Acceleration/Head Acceleraion)')
title('No Outliers')
hold off
suptitle('Constant Acceleration Gain AM Method')
savefig('ImpulseDataByGroupAM.fig')
saveas(gcf,'ImpulseDataByGroupAM.jpg')
cd ../

save('ImpulseFigureData.mat','-struct','rampfig')