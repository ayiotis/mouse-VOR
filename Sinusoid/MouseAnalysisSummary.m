%% Mouse Analysis Summary
%This function takes in a mouse ID and creates a table called mousetab with
%the Gain, phase, standard deviations for the parameters and total number
%of cycles that went into those numbers for each frequency. It also makes a
%figure that shows all of the data for that mouse and the sine fits for
%those data.

%Written by Andrianna Ayiotis
%Last updated 06/22/2018

function mousetab = MouseAnalysisSummary(mouse)
%This needs all the frequencies to have been analyzed already
freqs = {'10','5','2','1','0p5','0p2','0p1','0p05','0p02'};
freqs_num = str2double(strrep(freqs,'p','.'));
dat = zeros(length(freqs),5);
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(freqs)
    f_name = [mouse,'-',freqs{i},'HzSineAnalyzed.mat'];
    load(f_name,'SineAnalyzed')
    dat(i,:) = SineAnalyzed.sum_tab{3,:};
    %Make one big plot
    L = SineAnalyzed.L;
    R = SineAnalyzed.R;
    Chair = SineAnalyzed.Chair;
    t = SineAnalyzed.t; 
    maxvel = 1.2*max([max([Chair.cycle_avg,L.eyefit,R.eyefit]),abs(min([Chair.cycle_avg,L.eyefit,R.eyefit]))]);
    subplot(3,3,i)
    plot(t,L.All,'Color',[1,0.7,0.7]) %color is light red
    hold on   
    plot(t,R.All,'Color',[1,0.7,1]) %color is light pink
    p1 = plot(t,Chair.cycle_avg,'k');
    p3 = plot(t,L.Fit,'Color','r');
    p6 = plot(t,R.Fit,'Color','m'); 
    p4 = plot(t,L.eyefit,'r','LineWidth',2);
    p7 = plot(t,R.eyefit,'m','LineWidth',2);
    hold off
    title([strrep(freqs{i},'p','.'),' Hz'])
    axis([t(1),t(end),-maxvel,maxvel])
    if(i==7||i==8||i==9)
       xlabel('Time (s)') 
    end
    if(i==1||i==4||i==7)
       ylabel('Angular Velocity (dps)') 
    end
end
legend([p1,p3(1),p4,p6(1),p7],{'Inverted Chair','Left Traces','Left Fit','Right Traces','Right Fit'})
suptitle([mouse,' Sinusoidal Response'])
lab = horzcat({'Frequency'},SineAnalyzed.sum_tab.Properties.VariableNames);
mousetab = array2table([freqs_num',dat]);
mousetab.Properties.VariableNames = lab;
end