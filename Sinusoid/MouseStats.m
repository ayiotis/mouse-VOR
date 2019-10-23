%% Mouse Stats
%This script analyzes the data in the table in MouseSummary.mat

%Written by Andrianna Ayiotis
%Last updated 07/27/2018
%Last run on 07/27/2018 

%% Load the table
load('MouseSummary.mat','allmousetab')
freqs = allmousetab.Frequency;
freq = freqs(1:9);
cd Figures
%% Make Bode plots split by mouse type
het_tab = allmousetab(contains(allmousetab.Type,'het'),:);
cko_tab = allmousetab(contains(allmousetab.Type,'cko'),:);
het_mice = unique(het_tab.Mouse);
cko_mice = unique(cko_tab.Mouse);
het_gains = reshape(het_tab.Gain,9,[])';
het_phases = reshape(het_tab.Phase,9,[])';
cko_gains = reshape(cko_tab.Gain,9,[])';
cko_phases = reshape(cko_tab.Phase,9,[])';

figure
subplot(2,1,1)
loglog(freq,mean(het_gains),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none')
hold on
errorbar(freq,mean(het_gains),std(het_gains),'k')
loglog(freq,mean(cko_gains),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none')
errorbar(freq,mean(cko_gains),std(cko_gains),'b')
hold off
title('Frequency Response: Gain')
ylabel('Gain')
set(gca,'XLimMode','manual')
set(gca,'XLim',[0.01 11])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
set(gca,'YLimMode','manual')
set(gca,'YLim',[0.0001 10])
set(gca,'YTickMode','manual')
set(gca,'YTick',[0.0001 0.001 0.01 0.1 1 10])
set(gca,'YTickLabelMode','manual')
set(gca,'YTickLabel',num2cell([0.0001 0.001 0.01 0.1 1 10]))
subplot(2,1,2)
p3 = semilogx(freq,mean(het_phases),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');
hold on
errorbar(freq,mean(het_phases),std(het_phases),'k')
p4 = semilogx(freq,mean(cko_phases),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none');
errorbar(freq,mean(cko_phases),std(cko_phases),'b')
hold off
title('Frequency Response: Phase')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
legend([p3(1),p4(1)],{'Het Mice','Cko Mice'})
set(gca,'XLimMode','manual')
set(gca,'XLim',[0.01 11])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
savefig('GainPhaseDataByType.fig')
saveas(gcf,'GainPhaseDataByType.jpg')
%% Fit high pass filters to each mouse
%First order
H_eq1 = @(p,ff) 1j*p(2).*ff./(p(1)+1j.*ff);
%Second order single pole filter
H_eq2 = @(p,ff) -p(2).*(ff./(p(1)+1j.*ff)).^2;
%Fmincon params
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,0];
ub = [Inf,Inf];
het_params = zeros(4,length(het_mice));
for i = 1:length(het_mice)
    sub_het_tab = het_tab(contains(het_tab.Mouse,het_mice{i}),:);
    het_f = sub_het_tab.Frequency;
    het_g = sub_het_tab.Gain;
    het_g_sd = sub_het_tab.Gain_std;
    het_p = sub_het_tab.Phase;
    het_p_sd = sub_het_tab.Phase_std;
    het_n = sub_het_tab.n_cycles;    
    het_resp = het_g.*(cosd(het_p)+1j*sind(het_p));
    LSCF1 = @(p) sum(abs(H_eq1(p,het_f)-het_resp),'omitnan');
    het_params(1:2,i) = fmincon(LSCF1,[0.1;1],A,b,Aeq,beq,lb,ub);
    LSCF2 = @(p) sum(abs(H_eq2(p,het_f)-het_resp),'omitnan');
    het_params(3:4,i) = fmincon(LSCF2,[0.1;1],A,b,Aeq,beq,lb,ub);
end
het_params2 = het_params(:,2:end);
het_gains_rmo = het_gains(2:end,:);
het_phases_rmo = het_phases(2:end,:);
het_fit1_all = mean(het_params(1:2,:),2);
het_fit2_all = mean(het_params(3:4,:),2);
het_fit1_rmo = mean(het_params2(1:2,:),2);
het_fit2_rmo = mean(het_params2(3:4,:),2);
cko_params = zeros(4,length(cko_mice));
for i = 1:length(cko_mice)
    sub_cko_tab = cko_tab(contains(cko_tab.Mouse,cko_mice{i}),:);
    cko_f = sub_cko_tab.Frequency;
    cko_g = sub_cko_tab.Gain;
    cko_g_sd = sub_cko_tab.Gain_std;
    cko_p = sub_cko_tab.Phase;
    cko_p_sd = sub_cko_tab.Phase_std;
    cko_n = sub_cko_tab.n_cycles;    
    cko_resp = cko_g.*(cosd(cko_p)+1j*sind(cko_p));
    LSCF1 = @(p) sum(abs(H_eq1(p,cko_f)-cko_resp));
    cko_params(1:2,i) = fmincon(LSCF1,[0.1;1],A,b,Aeq,beq,lb,ub);
    LSCF2 = @(p) sum(abs(H_eq2(p,cko_f)-cko_resp));
    cko_params(3:4,i) = fmincon(LSCF2,[0.1;1],A,b,Aeq,beq,lb,ub);
end
cko_params2 = cko_params(:,[1:3,5]);
cko_gains_rmo = cko_gains([1:3,5],:);
cko_phases_rmo = cko_phases([1:3,5],:);
cko_fit1_all = mean(cko_params(1:2,:),2);
cko_fit2_all = mean(cko_params(3:4,:),2);
cko_fit1_rmo = mean(cko_params2(1:2,:),2);
cko_fit2_rmo = mean(cko_params2(3:4,:),2);
%% Summary Plots
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
loglog(freq,mean(het_gains),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none')
hold on
errorbar(freq,mean(het_gains),std(het_gains),'k.')
loglog(freq,abs(H_eq1(het_fit1_all,freq)),'k--')
loglog(freq,abs(H_eq2(het_fit2_all,freq)),'k-')
hold off
title('Frequency Response for All Data: Gain')
ylabel('Gain')
axis([0.01,11,0.00008,2])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
set(gca,'YTickMode','manual')
set(gca,'YTick',[0.001 0.01 0.1 1 10])
set(gca,'YTickLabelMode','manual')
set(gca,'YTickLabel',num2cell([0.001 0.01 0.1 1 10]))
subplot(2,2,3)
semilogx(freq,mean(het_phases),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none')
hold on
errorbar(freq,mean(het_phases),std(het_phases),'k.')
semilogx(freq,180/pi*angle(H_eq1(het_fit1_all,freq)),'k--')
semilogx(freq,180/pi*angle(H_eq2(het_fit2_all,freq)),'k-')
hold off
title('Frequency Response for All Data: Phase')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
axis([0.01,11,-30,175])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
subplot(2,2,2)
loglog(freq,mean(het_gains_rmo),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none')
hold on
errorbar(freq,mean(het_gains_rmo),std(het_gains_rmo),'k.')
loglog(freq,abs(H_eq1(het_fit1_rmo,freq)),'k--')
loglog(freq,abs(H_eq2(het_fit2_rmo,freq)),'k-')
hold off
title('Frequency Response with Outliers Removed: Gain')
axis([0.01,11,0.00008,2])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
set(gca,'YTickMode','manual')
set(gca,'YTick',[0.001 0.01 0.1 1])
set(gca,'YTickLabelMode','manual')
set(gca,'YTickLabel',num2cell([0.001 0.01 0.1 1]))
subplot(2,2,4)
p1 = semilogx(freq,mean(het_phases_rmo),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');
hold on
errorbar(freq,mean(het_phases_rmo),std(het_phases_rmo),'k.')
p2 = semilogx(freq,180/pi*angle(H_eq1(het_fit1_rmo,freq)),'k--');
p3 = semilogx(freq,180/pi*angle(H_eq2(het_fit2_rmo,freq)),'k-');
hold off
title('Frequency Response with Outliers Removed: Phase')
xlabel('Frequency (Hz)')
legend([p1(1),p2,p3],{'Mean Values','1st Order Fit','2nd Order Fit'});
axis([0.01,11,-30,175])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
suptitle('Het Mice Sinusoidal Responses')
savefig('HighPassFilterHetFitsAll.fig')
saveas(gcf,'HighPassFilterHetFitsAll.jpg')

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
loglog(freq,mean(cko_gains),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none')
hold on
errorbar(freq,mean(cko_gains),std(cko_gains),'b.')
loglog(freq,abs(H_eq1(cko_fit1_all,freq)),'b--')
loglog(freq,abs(H_eq2(cko_fit2_all,freq)),'b-')
hold off
title('Frequency Response for All Data: Gain')
ylabel('Gain')
axis([0.01 11 0.0009 2])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
set(gca,'YTickMode','manual')
set(gca,'YTick',[0.001 0.01 0.1 1])
set(gca,'YTickLabelMode','manual')
set(gca,'YTickLabel',num2cell([0.001 0.01 0.1 1]))
subplot(2,2,3)
semilogx(freq,mean(cko_phases),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none')
hold on
errorbar(freq,mean(cko_phases),std(cko_phases),'b.')
semilogx(freq,180/pi*angle(H_eq1(cko_fit1_all,freq)),'b--')
semilogx(freq,180/pi*angle(H_eq2(cko_fit2_all,freq)),'b-')
hold off
title('Frequency Response for All Data: Phase')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
axis([0.01,11,-35,165])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
subplot(2,2,2)
loglog(freq,mean(cko_gains_rmo),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none')
hold on
errorbar(freq,mean(cko_gains_rmo),std(cko_gains_rmo),'b.')
loglog(freq,abs(H_eq1(cko_fit1_rmo,freq)),'b--')
loglog(freq,abs(H_eq2(cko_fit2_rmo,freq)),'b-')
hold off
title('Frequency Response with Outliers Removed: Gain')
axis([0.01,11,0.0009,2])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
set(gca,'YTickMode','manual')
set(gca,'YTick',[0.001 0.01 0.1 1 10])
set(gca,'YTickLabelMode','manual')
set(gca,'YTickLabel',num2cell([0.001 0.01 0.1 1 10]))
subplot(2,2,4)
p1 = semilogx(freq,mean(cko_phases_rmo),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none');
hold on
errorbar(freq,mean(cko_phases_rmo),std(cko_phases_rmo),'b.')
p2 = semilogx(freq,180/pi*angle(H_eq1(cko_fit1_rmo,freq)),'b--');
p3 = semilogx(freq,180/pi*angle(H_eq2(cko_fit2_rmo,freq)),'b-');
hold off
title('Frequency Response with Outliers Removed: Phase')
xlabel('Frequency (Hz)')
legend([p1(1),p2,p3],{'Mean Values','1st Order Fit','2nd Order Fit'});
axis([0.01,11,-35,165])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
suptitle('Cko Mice Sinusoidal Response')
savefig('HighPassFilterCkoFitsAll.fig')
saveas(gcf,'HighPassFilterCkoFitsAll.jpg')

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
loglog(freq,mean(het_gains_rmo),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none')
hold on
errorbar(freq,mean(het_gains_rmo),std(het_gains_rmo),'k.')
loglog(freq,mean(cko_gains_rmo),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none')
errorbar(freq,mean(cko_gains_rmo),std(cko_gains_rmo),'b.')
loglog(freq,abs(H_eq1(het_fit1_rmo,freq)),'k-')
loglog(freq,abs(H_eq1(cko_fit1_rmo,freq)),'b-')
hold off
title('Frequency Response: Gain')
ylabel('Gain')
axis([0.01,11,0.00008,2])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
set(gca,'YTickMode','manual')
set(gca,'YTick',[0.0001 0.001 0.01 0.1 1])
set(gca,'YTickLabelMode','manual')
set(gca,'YTickLabel',num2cell([0.0001 0.001 0.01 0.1 1]))
subplot(2,1,2)
p1 = semilogx(freq,mean(het_phases_rmo),'Color','k','Marker','o','MarkerFaceColor','k','LineStyle','none');
hold on
errorbar(freq,mean(het_phases_rmo),std(het_phases_rmo),'k.');
p2 = semilogx(freq,mean(cko_phases_rmo),'Color','b','Marker','^','MarkerFaceColor','b','LineStyle','none');
errorbar(freq,mean(cko_phases_rmo),std(cko_phases_rmo),'b.');
p3 = semilogx(freq,180/pi*angle(H_eq1(het_fit1_rmo,freq)),'k-');
p4 = semilogx(freq,180/pi*angle(H_eq1(cko_fit1_rmo,freq)),'b-');
hold off
title('Frequency Response: Phase')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
legend([p1(1),p2(1),p3,p4],{'Het Avg','Cko Avg','Het 1st Order Fit','Cko 1st Order Fit'})
set(gca,'XLimMode','manual')
set(gca,'XLim',[0.01 11])
set(gca,'XTickMode','manual')
set(gca,'XTick',[0.01 0.1 1 10])
set(gca,'XTickLabelMode','manual')
set(gca,'XTickLabel',num2cell([0.01 0.1 1 10]))
savefig('HighPassFilterFits.fig')
saveas(gcf,'HighPassFilterFits.jpg')

sinfig.freq = freq;
sinfig.H_eq1 = H_eq1;
sinfig.H_eq2 = H_eq2;
sinfig.het_gains_rmo = het_gains_rmo;
sinfig.cko_gains_rmo = cko_gains_rmo;
sinfig.het_phases_rmo = het_phases_rmo;
sinfig.cko_phases_rmo = cko_phases_rmo;
sinfig.het_fit1_rmo = het_fit1_rmo;
sinfig.cko_fit1_rmo = cko_fit1_rmo;
sinfig.het_fit2_rmo = het_fit2_rmo;
sinfig.cko_fit2_rmo = cko_fit2_rmo;
%% All p-values
cd ../
save('SineFigureData','-struct','sinfig');
p1 = ranksum(het_params(1,:),cko_params(1,:)); %all fc 1st order
p2 = ranksum(het_params(2,:),cko_params(2,:)); %all K 1st order
p5 = ranksum(het_params2(1,:),cko_params2(1,:)); %no outliers fc 1st
p6 = ranksum(het_params2(2,:),cko_params2(2,:)); %no outliers K 1st
p9 = ranksum(het_params(3,:),cko_params(3,:)); %all fc 2nd
p10 = ranksum(het_params(4,:),cko_params(4,:)); %all K 2nd
p13 = ranksum(het_params2(3,:),cko_params2(3,:)); %no outliers fc 2nd
p14 = ranksum(het_params2(4,:),cko_params2(4,:)); %no outliers K 2nd
[~,p3] = ttest2(het_params(1,:),cko_params(1,:),'Vartype','unequal'); 
[~,p4] = ttest2(het_params(2,:),cko_params(2,:),'Vartype','unequal'); 
[~,p7] = ttest2(het_params2(1,:),cko_params2(1,:),'Vartype','unequal'); 
[~,p8] = ttest2(het_params2(2,:),cko_params2(2,:),'Vartype','unequal'); 
[~,p11] = ttest2(het_params(3,:),cko_params(3,:),'Vartype','unequal'); 
[~,p12] = ttest2(het_params(4,:),cko_params(4,:),'Vartype','unequal'); 
[~,p15] = ttest2(het_params2(3,:),cko_params2(3,:),'Vartype','unequal'); 
[~,p16] = ttest2(het_params2(4,:),cko_params2(4,:),'Vartype','unequal'); 
p = [p1,p2;p3,p4;p5,p6;p7,p8;p9,p10;p11,p12;p13,p14;p15,p16];
stat_sum = [cell2table({'Mann-Whitney U-test','All 1st Order';'Two-Sided t-test','All 1st Order';'Mann-Whitney U-test','No Outliers 1st Order';'Two-sided t-test','No Outliers 1st Order';'Mann-Whitney U-test','All 2nd Order';'Two-Sided t-test','All 2nd Order';'Mann-Whitney U-test','No Outliers 2nd Order';'Two-sided t-test','No Outliers 2nd Order'}),array2table(p)];
stat_sum.Properties.VariableNames = {'StatisticalTest','DataUsed','CenterFrequency','HighPassGain'};
disp('p-values:')
disp(stat_sum)
save('MouseStatisticsSummary.mat','stat_sum')