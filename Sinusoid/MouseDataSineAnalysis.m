%% Mouse Data Sine Analysis
%This script takes in the struct CleanDat made by MouseDataSineCleanData
%and calculates gain and phase with a programmatic sine fit or manual user
%selection.

%Written by Andrianna Ayiotis
%Last updated 08/06/2018

function SineAnalyzed = MouseDataSineAnalysis(CleanDat)
redo = 'y';
while(strcmp(redo,'y'))
    %% Load data from file
    round_freq = CleanDat.info.round_freq;
    mouse = CleanDat.info.mouse;
    t = CleanDat.t;
    m_chair = CleanDat.Chair.cycle_avg;
    Leye_seg = CleanDat.LEye.cycles;
    Reye_seg = CleanDat.REye.cycles;
    L = CleanDat.LEye;
    R = CleanDat.REye;
    Chair = CleanDat.Chair;

    maxvel = 150;
    chair_freq = CleanDat.Chair.freq;
    chair_amp = CleanDat.Chair.amp;
    %% Plot data
    p1 = plot(t,m_chair,'k','LineWidth',2);
    hold on
    p2 = plot(t,Leye_seg,'Color','r');
    p3 = plot(t,Reye_seg,'Color','m');
    hold off
    legend([p1,p2(1),p3(1)],{'Chair','Left Eye Traces','Right Eye Traces'})  
    title([mouse,' ',num2str(round_freq),' Hz Cycles'])
    xlabel('Time (s)')
    ylabel('Velocity (dps)')
    axis([t(1),t(end),-maxvel,maxvel])
    %% Find amplitude and phase
    disp('There are multiple methods you can use to analyze the sine waves. Select one to try.')
    disp('1. Whole Cycle Fit')
    disp('2. Half Cycle Fit')

    method = input('Choose the number you would like to use: ','s');
    while(~strcmp(method,'1')&&~strcmp(method,'2'))
       disp('You have entered an invalid selection. Enter only the numbers 1, or 2 with no quotation marks around them.')
       method = input('Choose the number you would like to use: ','s');
    end

    %Sine function used in all methods
    fit = @(p,tt)  p(1).*(sin(2*pi*chair_freq*tt + pi/180*p(2))) + p(3);
    options = optimset('MaxFunEvals',1000);

    if(strcmp(method,'1'))
        methods = 'Whole Cycle Least Squares Cost Function Minimization';
        if(~all(isnan(reshape(Leye_seg,[],1))))
            %Left eye
            L_params = zeros(3,size(Leye_seg,1));
            for i = 1:size(Leye_seg,1)
                %Make intitial guesses for amplitude, phase and offset (freq should be same
                %as the chair)
                L_ampg = (max(Leye_seg(i,:)) - min(Leye_seg(i,:)))/2;
                L_offg = max(Leye_seg(i,:)) - L_ampg;               
                %Definie the least squares cost function to minimize
                LSCF = @(p) sum((fit(p,t) - Leye_seg(i,:)).^2,'omitnan');   
                L_params(:,i) = fminsearch(LSCF, [L_ampg;0;L_offg],options);
                %Make sure amplitude isn't negative and update phase if it
                %is
                if(L_params(1,i)<0)
                   L_params(1,i) = abs(L_params(1,i));
                   L_params(2,i) = L_params(2,i)+180;
                end            
            end
            %Now fit the whole trace
            long_t = repmat(t,1,size(Leye_seg,1));
            LSCF = @(p) sum((fit(p,long_t) - reshape(Leye_seg',1,[])).^2,'omitnan');   
            L_params_avg = fminsearch(LSCF,[max(reshape(Leye_seg,1,[]));0;0],options);
            if(L_params_avg(1)<0)
                L_params_avg(1) = abs(L_params_avg(1));
                L_params_avg(2) = L_params_avg(2)+180;
            end
            L_params_std = std(L_params,[],2,'omitnan');
            L_eyefit = fit(L_params_avg,t);
            %Make a summary of all the relevant values
            L_vals = [L_params_avg(1)/chair_amp,L_params_std(1)/chair_amp,L_params_avg(2),L_params_std(2),size(Leye_seg,1)];
            L.eyefit = L_eyefit;
            L.params = L_params; %amplitude, phase, vertical shift
            L.params_avg = L_params_avg;
            L.params_std = L_params_std;
        else
            L_vals = [nan,nan,nan,nan,0];
            L.eyefit = NaN(1,length(t));
            L.params = [nan;nan;nan]; %amplitude, phase, vertical shift
            L.params_avg = [nan;nan;nan];
            L.params_std = [nan;nan;nan];
        end
        if(~all(isnan(reshape(Reye_seg,[],1))))
            %Right eye
            R_params = zeros(3,size(Reye_seg,1));
            for i = 1:size(Reye_seg,1)
                %Make intitial guesses for amplitude, phase and offset (freq should be same
                %as the chair)
                R_ampg = (max(Reye_seg(i,:)) - min(Reye_seg(i,:)))/2;
                R_offg = max(Reye_seg(i,:)) - R_ampg;               
                %Definie the least squares cost function to minimize
                RSCF = @(p) sum((fit(p,t) - Reye_seg(i,:)).^2,'omitnan');   
                R_params(:,i) = fminsearch(RSCF, [R_ampg;0;R_offg],options);
                %Make sure amplitude isn't negative
                if(R_params(1,i)<0)
                   R_params(1,i) = abs(R_params(1,i));
                   R_params(2,i) = R_params(2,i)+180;
                end               
            end
            %Now fit the whole trace
            long_t = repmat(t,1,size(Reye_seg,1));
            LSCF = @(p) sum((fit(p,long_t) - reshape(Reye_seg',1,[])).^2,'omitnan');   
            R_params_avg = fminsearch(LSCF,[max(reshape(Reye_seg,1,[]));0;0],options);
            if(R_params_avg(1)<0)
                R_params_avg(1) = abs(R_params_avg(1));
                R_params_avg(2) = R_params_avg(2)+180;
            end
            R_params_std = std(R_params,[],2,'omitnan');
            R_eyefit = fit(R_params_avg,t);
            %Summary vector
            R_vals = [R_params_avg(1)/chair_amp,R_params_std(1)/chair_amp,R_params_avg(2),R_params_std(2),size(Reye_seg,1)];
            R.eyefit = R_eyefit;
            R.params = R_params; %amplitude, phase, vertical shift
            R.params_avg = R_params_avg;
            R.params_std = R_params_std;
        else 
            R_vals = [nan,nan,nan,nan,0];
            R.eyefit = NaN(1,length(t));
            R.params = [nan;nan;nan]; %amplitude, phase, vertical shift
            R.params_avg = [nan;nan;nan];
            R.params_std = [nan;nan;nan];
        end
    elseif(strcmp(method,'2'))
        methods = 'Half Cycle Fit Least Squares Cost';
        if(~all(isnan(reshape(Leye_seg,[],1))))
            %Left eye
            %Plot to choose which side
            plot(t,m_chair,'k','LineWidth',2);
            hold on
            plot(t,Leye_seg,'r'); 
            hold off
            side_input = input('Which half cycle would you like to use (positive=p, negative=n)? (p/n) ','s');
            while(~strcmp(side_input,'p')&&~strcmp(side_input,'n'))
               disp('Invalid character entered. Enter only "p" or "n" without the quotation marks.')
               side_input = input('Which half cycle would you like to use (positive=p, negative=n)? (p/n) ','s');
            end
            mid = floor(length(t)/2);
            if(strcmp(side_input,'p'))
                Leye_seg(:,mid:end) = nan;
                methods = [methods,' Left Eye Positive Cycle'];
            else
                Leye_seg(:,1:mid) = nan;
                methods = [methods,' Left Eye Positive Cycle'];
            end          
            L_params = zeros(3,size(Leye_seg,1));   
            for i = 1:size(Leye_seg,1)
                %Make intitial guesses for amplitude, phase and offset (freq should be same
                %as the chair)
                L_ampg = abs(max(Leye_seg(i,:),[],'omitnan') - min(Leye_seg(i,:),[],'omitnan'));
                %Definie the least squares cost function to minimize
                LSCF = @(p) sum((fit(p,t) - Leye_seg(i,:)).^2,'omitnan');                             
                L_params(:,i) = fminsearch(LSCF, [L_ampg;0;0],options);
                if(L_params(1,i)<0)
                   L_params(1,i) = abs(L_params(1,i));
                   L_params(2,i) = L_params(2,i)+180;
                end              
            end
            %Now fit the whole trace
            long_t = repmat(t,1,size(Leye_seg,1));
            LSCF = @(p) sum((fit(p,long_t) - reshape(Leye_seg',1,[])).^2,'omitnan');   
            L_params_avg = fminsearch(LSCF,[max(reshape(Leye_seg,1,[]));0;0],options);
            if(L_params_avg(1)<0)
                L_params_avg(1) = abs(L_params_avg(1));
                L_params_avg(2) = L_params_avg(2)+180;
            end
            L_params_std = std(L_params,[],2,'omitnan');         
            L_eyefit = fit(L_params_avg,t);
            %Make a summary of all the relevant values
            L_vals = [L_params_avg(1)/chair_amp,L_params_std(1)/chair_amp,L_params_avg(2),L_params_std(2),size(Leye_seg,1)];
            L.eyefit = L_eyefit;
            L.params = L_params; %amplitude, phase, vertical shift
            L.params_avg = L_params_avg;
            L.params_std = L_params_std;
        else
            L_vals = [nan,nan,nan,nan,0];
            L.eyefit = NaN(1,length(t));
            L.params = [nan;nan;nan]; %amplitude, phase, vertical shift
            L.params_avg = [nan;nan;nan];
            L.params_std = [nan;nan;nan];
        end
        if(~all(isnan(reshape(Reye_seg,[],1))))
            %Right eye
            %Plot to choose side
            plot(t,m_chair,'k','LineWidth',2);
            hold on
            plot(t,Reye_seg,'Color','m'); 
            hold off
            side_input = input('Which half cycle would you like to use (positive=p, negative=n)? (p/n) ','s');
            while(~strcmp(side_input,'p')&&~strcmp(side_input,'n'))
               disp('Invalid character entered. Enter only "p" or "n" without the quotation marks.')
               side_input = input('Which half cycle would you like to use (positive=p, negative=n)? (p/n) ','s');
            end
            mid = floor(length(t)/2);
            if(strcmp(side_input,'p'))
                Reye_seg(:,mid:end) = nan;
                methods = [methods,' Right Eye Positive Cycle'];
            else
                Reye_seg(:,1:mid) = nan;
                methods = [methods,' Right Eye Negative Cycle'];
            end
            R_params = zeros(3,size(Reye_seg,1));           
            for i = 1:size(Reye_seg,1)
                %Make intitial guess for amplitude
                R_ampg = abs(max(Reye_seg(i,:),[],'omitnan') - min(Reye_seg(i,:),[],'omitnan'));
                
                %Definie the least squares cost function to minimize
                RSCF = @(p) sum((fit(p,t) - Reye_seg(i,:)).^2,'omitnan');                             
                R_params(:,i) = fminsearch(RSCF, [R_ampg;0;0],options);
                if(R_params(1,i)<0)
                   R_params(1,i) = abs(R_params(1,i));
                   R_params(2,i) = R_params(2,i)+180;
                end              
            end
            %Now fit the whole trace
            long_t = repmat(t,1,size(Reye_seg,1));
            LSCF = @(p) sum((fit(p,long_t) - reshape(Reye_seg',1,[])).^2,'omitnan');   
            R_params_avg = fminsearch(LSCF,[max(reshape(Reye_seg,1,[]));0;0],options);
            if(R_params_avg(1)<0)
                R_params_avg(1) = abs(R_params_avg(1));
                R_params_avg(2) = R_params_avg(2)+180;
            end
            R_params_std = std(R_params,[],2,'omitnan');
            R_eyefit = fit(R_params_avg,t);
            %Summary vector
            R_vals = [R_params_avg(1)/chair_amp,R_params_std(1)/chair_amp,R_params_avg(2),R_params_std(2),size(Reye_seg,1)];
            R.eyefit = R_eyefit;
            R.params = R_params; %amplitude, phase, vertical shift
            R.params_avg = R_params_avg;
            R.params_std = R_params_std;
        else 
            R_vals = [nan,nan,nan,nan,0];
            R.eyefit = NaN(1,length(t));
            R.params = [nan;nan;nan]; %amplitude, phase, vertical shift
            R.params_avg = [nan;nan;nan];
            R.params_std = [nan;nan;nan];
        end
    end
    %Combine the eyes
    if(all(isnan(reshape(Leye_seg,[],1)))&&all(isnan(reshape(Reye_seg,[],1)))) %If they're all nan
        Combo = L_vals; %This will also be NaN
    else 
        pool_gain_std = sqrt(sum([L_vals(2)^2*(L_vals(5)-1),R_vals(2)^2*(R_vals(5)-1)],'omitnan')/(L_vals(5)+R_vals(5)-2));
        pool_phase_std = sqrt(sum([L_vals(4)^2*(L_vals(5)-1),R_vals(4)^2*(R_vals(5)-1)],'omitnan')/(L_vals(5)+R_vals(5)-2));
        Combo = [mean([L_vals(1),R_vals(1)],'omitnan'),pool_gain_std,mean([L_vals(3),R_vals(3)],'omitnan'),pool_phase_std,L_vals(5)+R_vals(5)];
    end
    L.Fit = Leye_seg;
    R.Fit = Reye_seg;
    % Summary Plot    
    p2 = plot(t,L.All,'Color',[1,0.7,0.7]); %color is light red
    hold on 
    p5 = plot(t,R.All,'Color',[1,0.7,1]); %color is light pink
    p1 = plot(t,Chair.cycle_avg,'k');
    p3 = plot(t,L.Fit,'Color','r');
    p6 = plot(t,R.Fit,'Color','m'); 
    p4 = plot(t,L.eyefit,'r','LineWidth',8);
    p7 = plot(t,R.eyefit,'m','LineWidth',8);
    legend([p1,p2(1),p3(1),p4,p5(1),p6(1),p7],{'Chair','Raw Left Traces','Analyzed Left Traces','Left Eye Fit','Raw Right Traces','Analyzed Right Traces','Right Eye Fit'})
    hold off
    title([mouse,' ',num2str(round_freq),' Hz'])
    axis([t(1),t(end),-maxvel,maxvel])

    %Let the user redo if they would like to analyze the data with a
    %different method
    redo = input('Would you like to redo the analysis? (y/n) ','s');
    while(~strcmp(redo,'y')&&~strcmp(redo,'n'))
       disp('Invalid character entered. Enter only "y" or "n" without the quotation marks.')
       redo = input('Would you like to redo the analysis? (y/n) ','s');
    end 
end
%% Make summary tables and structs to store all necessary parameters
sum_tab = array2table([L_vals;R_vals;Combo]);
sum_tab.Properties.VariableNames = {'Gain','Gain_std','Phase','Phase_std','n_cycles'};
sum_tab.Properties.RowNames = {'LeftEye','RightEye','Combined'};
%Update and add to the info section
info = CleanDat.info;
info.CleanDat = CleanDat.info;
info.CleanDat = rmfield(info.CleanDat,'SegDat');
info = rmfield(info,'rm_qp');
info = rmfield(info,'qp_a_thresh_l');
info = rmfield(info,'qp_a_thresh_r');
info = rmfield(info,'maxvel');
info.infile = info.fname;
info.analysis_time = clock;
fname = [info.fname(1:end-18),'Analyzed.mat'];
info.fname = fname;
info.analysis_method = methods; 
SineAnalyzed.info = info;
SineAnalyzed.sum_tab = sum_tab;
SineAnalyzed.t = t;
SineAnalyzed.Chair = Chair;
SineAnalyzed.L = L;
SineAnalyzed.R = R;
end