%% Mouse Data Impulse Analysis
% This fucntion takes in the output of MouseDataImpulseCleanData, a struct 
% called CleanImpulseData that contains the traces fit for analysis and the
% parameters used to segment the code thus far.

% This script isolates the regions of interest for a linear fit of both
% chair and eye velocities and from those fits calculates gains and
% latencies for each trace. When plot_check = true, the user can chose to 
% see each individual fit and adjust the parameters if needed.

%Made by Andrianna Ayiotis
%Last updated 07/02/2018
function ImpulseAnalysisOutput = MouseDataImpulseAnalysis(CleanImpulseData,plot_check)
%% User inputs 
%Control what percents of the max/min is in the linear range for each eye.
lpLl = 0.3;
upLl = 0.8;
lpLr = 0.3;
upLr = 0.8;
lpRl = 0.3;
upRl = 0.8;
lpRr = 0.3;
upRr = 0.8;
%% Initialize 
t = CleanImpulseData.t;
Ll = CleanImpulseData.Ll;
Lr = CleanImpulseData.Lr;
Rl = CleanImpulseData.Rl;
Rr = CleanImpulseData.Rr;
ChairLl = Ll.chair_clean;
ChairRl = Rl.chair_clean;
ChairLr = Lr.chair_clean;
ChairRr = Rr.chair_clean;
EyeLl = Ll.eye_clean;
EyeRl = Rl.eye_clean;
EyeLr = Lr.eye_clean;
EyeRr = Rr.eye_clean;
info = CleanImpulseData.info;
minvel = info.minvel;
maxvel = info.maxvel;
mouse = info.mouse;
ms_p_ind = info.ms_p_ind;
targvel = maxvel - minvel;
subplot(1,1,1)
suptitle('')
%% Analyze the data to find gains/latencies for each of the curves (left eye going left)
if(~all(isnan(ChairLl)))
    %Chair fit on mean trace
    chair1 = mean(ChairLl);
    i1 = find(chair1>minvel&chair1<targvel);
    i2 = find(diff(i1)>1);
    t_chair = t(i1(1):i1(i2));
    chair2 = chair1(i1(1):i1(i2))';
    t_chair2 = [ones(length(chair2),1),t_chair];
    b = t_chair2\chair2;
    t_chair3 = [ones(length(t),1),t];
    chairfit = t_chair3*b;
    %Keep running this block of code until user confirmation
    u_rang = 'n';
    l_rang = 'n';
    while(strcmp(u_rang,'n')||strcmp(l_rang,'n')) 
        %Eye fit trace by trace
        GaLl = zeros(size(ChairLl,1),1);
        latLl = zeros(size(ChairLl,1),1); 
        AMGaLl = zeros(size(ChairLl,1),1);
        for i = 1:size(EyeLl,1)
            Lz1 = EyeLl(i,:);  
            i3 = find(chair1>1&chair1<=maxvel);
            i4 = find(diff(i3)>1);
            t_eye1 = t(i3(1)-1:i3(i4(1)+1)+floor(100/ms_p_ind));
            Lz2 = Lz1(i3(1)-1:i3(i4(1)+1)+floor(100/ms_p_ind));
            i5 = find(Lz2>lpLl*max(Lz2) & Lz2<upLl*max(Lz2));
            i6 = find(diff(i5)>1);
            if(isempty(i6))
                t_eye2 = t_eye1(i5(1):i5(end));
                Lz3 = Lz2(i5(1):i5(end))';
            else
                t_eye2 = t_eye1(i5(1):i5(i6));
                Lz3 = Lz2(i5(1):i5(i6))';
            end
            i7 = find(diff(Lz3)<0);
            if(~isempty(i7))
                t_eye1 = t(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                Lz2 = Lz1(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                i5 = find(Lz2>lpLl*max(Lz2) & Lz2<upLl*max(Lz2));
                i6 = find(diff(i5)>1);  
                if(isempty(i6))
                    t_eye2 = t_eye1(i5(1):i5(end));
                    Lz3 = Lz2(i5(1):i5(end))';
                else
                    t_eye2 = t_eye1(i5(1):i5(i6));
                    Lz3 = Lz2(i5(1):i5(i6))';
                end
            end
            t_eye3 = [ones(length(t_eye2),1),t_eye2];
            b2 = t_eye3\Lz3;
            eyefit = t_chair3*b2;
            %Gain and latency calculations
            GaLl(i,1) = b2(2)/b(2); %ratio
            latLl(i,1) = 1000*(-b2(1)/b2(2)+b(1)/b(2)); %in ms        
            %Calculate gain the AM way without a latency shift
            i_s = find(t==t_eye2(1));
            AMGaLl(i,1) = mean(eyefit(i_s:i_s+length(Lz3)-1)./chairfit(i_s:i_s+length(Lz3)-1));      
            %Plot
            if(plot_check)
                plot(t,chair1,'k',t,Lz1,'r',t,chairfit,'k',t,eyefit,'r',t(i_s:i_s+length(Lz3)-1),eyefit(i_s:i_s+length(Lz3)-1),'r*')
                title(['Check Fits for Left Eye Going to the Left on ',mouse])
                xlabel('Time (s)')
                ylabel('Angular Velocity (dps)')
                legend('Inverted Chair','Left Eye','Chair Fit','Left Eye Fit')
                axis([t(1) t(end) -0.5*maxvel 1.3*maxvel])
                pause;
            end
        end
        m_GaLl = mean(GaLl);
        m_latLl = mean(latLl);
        b3 = [m_GaLl*b(2)*(b(1)/b(2)-m_latLl/1000);m_GaLl*b(2)];
        m_eyefit = t_chair3*b3;
        % Show summary plot
        p1 = plot(t,chair1,'k');
        hold on
        p2 = plot(t,EyeLl,'r');
        p3 = plot(t,chairfit,'k','LineWidth',2);
        p4 = plot(t,m_eyefit,'r','LineWidth',2);
        title(['Check Fits for Left Eye Going to the Left on ',mouse])
        xlabel('Time (s)')
        ylabel('Angular Velocity (dps)')
        legend([p1,p2(1),p3,p4],{'Inverted Chair','Left Eye','Chair Fit','Left Eye Fit'})
        axis([t(1) t(end) -0.5*maxvel 1.3*maxvel])
        hold off
        %Check if need to rerun at all
        rerun = input('Is the fit appropriate? (y/n): ','s');
        while ~(strcmp(rerun,'y')||strcmp(rerun,'n')) 
            disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
            rerun = input('Is the fit appropriate? (y/n): ','s');
        end
        if(strcmp(rerun,'n'))
            %Check the upper linear level
            u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            while ~(strcmp(u_rang,'y')||strcmp(u_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            end
            if(strcmp(u_rang,'n'))
                u_str = ['The current upper limit percentage is ',num2str(upLl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upLl = new_rang/100;
            end
            %Check the lower linear level
            l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            while ~(strcmp(l_rang,'y')||strcmp(l_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            end
            if(strcmp(l_rang,'n'))
                u_str = ['The current lower limit percentage is ',num2str(lpLl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpLl = new_rang/100;
            end
            %Make sure the limits are valid
            while(lpLl>=upLl)
                %Reset upper limit
                disp('The lower percentage limit cannot be greater than or equal to the upper percentage limit. They must be changed.')
                u_str = ['The current upper limit percentage is ',num2str(upLl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upLl = new_rang/100;
                %Reset lower limit
                u_str = ['The current lower limit percentage is ',num2str(lpLl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpLl = new_rang/100;
            end
        else
            u_rang = 'y';
            l_rang = 'y';
        end
    end
    Ll.chair_m = chair1;
    Ll.chairfit = chairfit;
    Ll.eyefit = m_eyefit;
else
    disp('No traces accepted for left eye going left')
    GaLl = [];
    latLl = [];
    AMGaLl = [];
    Ll.chair_m = NaN(1,length(t));
    Ll.chairfit = NaN(1,length(t));
    Ll.eyefit = NaN(1,length(t));
end
%% Analyze the data to find gains/latencies for each of the curves (left eye going right)
if(~all(isnan(ChairLr)))
    %Chair fit on mean trace
    chair1 = mean(ChairLr);
    i1 = find(chair1<-minvel & chair1>-targvel);
    i2 = find(diff(i1)>1);
    t_chair = t(i1(1):i1(i2));
    chair2 = chair1(i1(1):i1(i2))';
    t_chair2 = [ones(length(chair2),1),t_chair];
    b = t_chair2\chair2;
    t_chair3 = [ones(length(t),1),t];
    chairfit = t_chair3*b;
    %Keep running this block of code until user confirmation
    u_rang = 'n';
    l_rang = 'n';
    while(strcmp(u_rang,'n')||strcmp(l_rang,'n')) 
        %Eye fit trace by trace
        GaLr = zeros(size(ChairLr,1),1);
        latLr = zeros(size(ChairLr,1),1);
        AMGaLr = zeros(size(ChairLr,1),1);
        for i = 1:size(ChairLr,1)
            Lz1 = EyeLr(i,:);
            i3 = find(chair1<-1&chair1>=-maxvel);
            i4 = find(diff(i3)>1);
            t_eye1 = t(i3(1)-1:i3(i4(1)+1)+floor(50/ms_p_ind));
            Lz2 = Lz1(i3(1)-1:i3(i4(1)+1)+floor(50/ms_p_ind));
            i5 = find(Lz2<lpLr*min(Lz2) & Lz2>upLr*min(Lz2));
            i6 = find(diff(i5)>1);
            if(isempty(i6))
                t_eye2 = t_eye1(i5(1):i5(end));
                Lz3 = Lz2(i5(1):i5(end))';
            else
                t_eye2 = t_eye1(i5(1):i5(i6));
                Lz3 = Lz2(i5(1):i5(i6))';
            end
            i7 = find(diff(Lz3)>0);
            if(~isempty(i7))
                t_eye1 = t(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                Lz2 = Lz1(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                i5 = find(Lz2<lpLr*min(Lz2) & Lz2>upLr*min(Lz2));
                i6 = find(diff(i5)>1);
                if(isempty(i6))
                    t_eye2 = t_eye1(i5(1):i5(end));
                    Lz3 = Lz2(i5(1):i5(end))';
                else
                    t_eye2 = t_eye1(i5(1):i5(i6));
                    Lz3 = Lz2(i5(1):i5(i6))';
                end
            end
            t_eye3 = [ones(length(t_eye2),1),t_eye2];
            b2 = t_eye3\Lz3;
            eyefit = t_chair3*b2;
            %Gain and latency calculations
            GaLr(i,1) = b2(2)/b(2); %ratio
            latLr(i,1) = 1000*(-b2(1)/b2(2) + b(1)/b(2)); %in ms
            %Calculate gain the AM way without latency shift
            i_s = find(t==t_eye2(1));
            AMGaLr(i,1) = mean(eyefit(i_s:i_s+length(Lz3)-1)./chairfit(i_s:i_s+length(Lz3)-1));
            %Plot
            if(plot_check)
                plot(t,chair1,'k',t,Lz1,'r',t,chairfit,'k',t,eyefit,'r',t(i_s:i_s+length(Lz3)-1),eyefit(i_s:i_s+length(Lz3)-1),'r*')
                title(['Check Fits for Left Eye Going to the Right on ',mouse])
                xlabel('Time (s)')
                ylabel('Angular Velocity (dps)')
                legend('Inverted Chair','Left Eye','Chair Fit','Left Eye Fit')
                axis([t(1) t(end) -1.3*maxvel 0.5*maxvel])
                pause;
            end
        end
        m_GaLr = mean(GaLr);
        m_latLr = mean(latLr);
        b3 = [m_GaLr*b(2)*(b(1)/b(2)-m_latLr/1000);m_GaLr*b(2)];
        m_eyefit = t_chair3*b3;
        %Show summary plots
        p1 = plot(t,chair1,'k');
        hold on
        p2 = plot(t,EyeLr,'r');
        p3 = plot(t,chairfit,'k','LineWidth',2);
        p4 = plot(t,m_eyefit,'r','LineWidth',2);
        title(['Check Fits for Left Eye Going to the Right on ',mouse])
        xlabel('Time (s)')
        ylabel('Angular Velocity (dps)')
        legend([p1,p2(1),p3,p4],{'Inverted Chair','Left Eye','Chair Fit','Left Eye Fit'})
        axis([t(1) t(end) -1.3*maxvel 0.5*maxvel])
        hold off
        %Check if need to rerun at all
        rerun = input('Is the fit appropriate? (y/n): ','s');
        while ~(strcmp(rerun,'y')||strcmp(rerun,'n')) 
            disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
            rerun = input('Is the fit appropriate? (y/n): ','s');
        end
        if(strcmp(rerun,'n'))
            %Check the upper linear level
            u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            while ~(strcmp(u_rang,'y')||strcmp(u_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            end
            if(strcmp(u_rang,'n'))
                u_str = ['The current upper limit percentage is ',num2str(upLr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upLr = new_rang/100;
            end
            %Check the lower linear level
            l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            while ~(strcmp(l_rang,'y')||strcmp(l_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            end
            if(strcmp(l_rang,'n'))
                u_str = ['The current lower limit percentage is ',num2str(lpLr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpLr = new_rang/100;
            end
            %Make sure the limits are valid
            while(lpLr>=upLr)
                disp('The lower percentage limit cannot be greater than or equal to the upper percentage limit. They must be changed.')
                u_str = ['The current upper limit percentage is ',num2str(upLr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upLr = new_rang/100;
                u_str = ['The current lower limit percentage is ',num2str(lpLr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpLr = new_rang/100;
            end
        else
            u_rang = 'y';
            l_rang = 'y';
        end
    end
    Lr.chair_m = chair1;
    Lr.chairfit = chairfit;
    Lr.eyefit = m_eyefit;
else
    disp('No traces accepted for left eye going right.')
    GaLr = [];
    latLr = [];
    AMGaLr = [];
    Lr.chair_m = NaN(1,length(t));
    Lr.chairfit = NaN(1,length(t));
    Lr.eyefit = NaN(1,length(t));
end
%% Analyze the data to find gains/latencies for each of the curves (right eye going left)
if(~all(isnan(ChairRl)))
    %Chair fit on mean trace
    chair1 = mean(ChairRl);
    i1 = find(chair1>minvel & chair1<targvel);
    i2 = find(diff(i1)>1);
    t_chair = t(i1(1):i1(i2));
    chair2 = chair1(i1(1):i1(i2))';
    t_chair2 = [ones(length(chair2),1),t_chair];
    b = t_chair2\chair2;
    t_chair3 = [ones(length(t),1),t];
    chairfit = t_chair3*b;
    %Keep running this block of code until user confirmation
    u_rang = 'n';
    l_rang = 'n';
    while(strcmp(u_rang,'n')||strcmp(l_rang,'n')) 
        %Eye fit trace by trace
        GaRl = zeros(size(ChairRl,1),1);
        latRl = zeros(size(ChairRl,1),1);
        AMGaRl = zeros(size(ChairRl,1),1);
        for i = 1:size(ChairRl,1)
            Rz1 = EyeRl(i,:);
            i3 = find(chair1>1&chair1<=maxvel);
            i4 = find(diff(i3)>1);
            t_eye1 = t(i3(1)-1:i3(i4(1)+1)+floor(50/ms_p_ind));
            Rz2 = Rz1(i3(1)-1:i3(i4(1)+1)+floor(50/ms_p_ind));
            i5 = find(Rz2>lpRl*max(Rz2) & Rz2<upRl*max(Rz2));
            i6 = find(diff(i5)>1);
            if(isempty(i6))
                t_eye2 = t_eye1(i5(1):i5(end));
                Rz3 = Rz2(i5(1):i5(end))';
            else
                t_eye2 = t_eye1(i5(1):i5(i6));
                Rz3 = Rz2(i5(1):i5(i6))';
            end
            i7 = find(diff(Rz3)<0);
            if(~isempty(i7))
                t_eye1 = t(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                Rz2 = Rz1(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                i5 = find(Rz2>lpRl*max(Rz2) & Rz2<upRl*max(Rz2));
                i6 = find(diff(i5)>1);
                if(isempty(i6))
                    t_eye2 = t_eye1(i5(1):i5(end));
                    Rz3 = Rz2(i5(1):i5(end))';
                else
                    t_eye2 = t_eye1(i5(1):i5(i6));
                    Rz3 = Rz2(i5(1):i5(i6))';
                end
            end
            t_eye3 = [ones(length(t_eye2),1),t_eye2];
            b2 = t_eye3\Rz3;
            eyefit = t_chair3*b2;
            %Gain and latency calculations
            GaRl(i,1) = b2(2)/b(2); %ratio
            latRl(i,1) = 1000*(-b2(1)/b2(2) + b(1)/b(2)); %in ms
            %Calculate gain the AM way wihtout latency shift
            i_s = find(t==t_eye2(1));
            AMGaRl(i,1) = mean(eyefit(i_s:i_s+length(Rz3)-1)./chairfit(i_s:i_s+length(Rz3)-1));
            %Plot
            if(plot_check)
                plot(t,chair1,'k',t,Rz1,'m',t,chairfit,'k',t,eyefit,'m',t(i_s:i_s+length(Rz3)-1),eyefit(i_s:i_s+length(Rz3)-1),'m*')
                title(['Check Fits for Right Eye Going to the Left on ',mouse])
                xlabel('Time (s)')
                ylabel('Angular Velocity (dps)')
                legend('Inverted Chair','Right Eye','Chair Fit','Right Eye Fit')
                axis([t(1) t(end) -0.5*maxvel 1.3*maxvel])
                pause;
            end
        end
        m_GaRl = mean(GaRl);
        m_latRl = mean(latRl);
        b3 = [m_GaRl*b(2)*(b(1)/b(2)-m_latRl/1000);m_GaRl*b(2)];
        m_eyefit = t_chair3*b3;
        p1 = plot(t,chair1,'k');
        hold on
        p2 = plot(t,EyeRl,'m');
        p3 = plot(t,chairfit,'k','LineWidth',2);
        p4 = plot(t,m_eyefit,'m','LineWidth',2);
        title(['Check Fits for Right Eye Going to the Left on ',mouse])
        xlabel('Time (s)')
        ylabel('Angular Velocity (dps)')
        legend([p1,p2(1),p3,p4],{'Inverted Chair','Right Eye','Chair Fit','Right Eye Fit'})
        axis([t(1) t(end) -0.5*maxvel 1.3*maxvel])
        hold off
        %Check if need to rerun at all
        rerun = input('Is the fit appropriate? (y/n): ','s');
        while ~(strcmp(rerun,'y')||strcmp(rerun,'n')) 
            disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
            rerun = input('Is the fit appropriate? (y/n): ','s');
        end
        if(strcmp(rerun,'n'))
            %Check the upper linear level
            u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            while ~(strcmp(u_rang,'y')||strcmp(u_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            end
            if(strcmp(u_rang,'n'))
                u_str = ['The current upper limit percentage is ',num2str(upRl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upRl = new_rang/100;
            end
            %Check the lower linear level
            l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            while ~(strcmp(l_rang,'y')||strcmp(l_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            end
            if(strcmp(l_rang,'n'))
                u_str = ['The current lower limit percentage is ',num2str(lpRl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpRl = new_rang/100;
            end
            %Make sure the limits are valid
            while(lpRl>=upRl)
                disp('The lower percentage limit cannot be greater than or equal to the upper percentage limit. They must be changed.')
                u_str = ['The current upper limit percentage is ',num2str(upRl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upRl = new_rang/100;
                u_str = ['The current lower limit percentage is ',num2str(lpRl*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpRl = new_rang/100;
            end
        else
            u_rang = 'y';
            l_rang = 'y';
        end
    end
    Rl.chair_m = chair1;
    Rl.chairfit = chairfit;
    Rl.eyefit = m_eyefit;
else
    disp('No traces accepted for right eye going left.')
    GaRl = [];
    latRl = [];
    AMGaRl = [];
    Rl.chair_m = NaN(1,length(t));
    Rl.chairfit = NaN(1,length(t));
    Rl.eyefit = NaN(1,length(t));
end
%% Analyze the data to find gains/latencies for each of the curves (right eye going right)
if(~all(isnan(ChairRr)))
    %Chair fit on mean trace
    chair1 = mean(ChairRr);
    i1 = find(chair1<-minvel & chair1>-targvel);
    i2 = find(diff(i1)>1);
    t_chair = t(i1(1):i1(i2));
    chair2 = chair1(i1(1):i1(i2))';
    t_chair2 = [ones(length(chair2),1),t_chair];
    b = t_chair2\chair2;
    t_chair3 = [ones(length(t),1),t];
    chairfit = t_chair3*b;
    %Keep running this block of code until user confirmation
    u_rang = 'n';
    l_rang = 'n';
    while(strcmp(u_rang,'n')||strcmp(l_rang,'n')) 
        %Eye fit trace by trace
        GaRr = zeros(size(ChairRr,1),1);
        latRr = zeros(size(ChairRr,1),1);
        AMGaRr = zeros(size(ChairRr,1),1);
        for i = 1:size(EyeRr,1)
            Rz1 = EyeRr(i,:);
            i3 = find(chair1<-1&chair1>=-maxvel);
            i4 = find(diff(i3)>1);
            t_eye1 = t(i3(1)-1:i3(i4(1)+1)+floor(50/ms_p_ind));
            Rz2 = Rz1(i3(1)-1:i3(i4(1)+1)+floor(50/ms_p_ind));
            i5 = find(Rz2<lpRr*min(Rz2) & Rz2>upRr*min(Rz2));
            i6 = find(diff(i5)>1);
            if(isempty(i6))
                t_eye2 = t_eye1(i5(1):i5(end));
                Rz3 = Rz2(i5(1):i5(end))';
            else
                t_eye2 = t_eye1(i5(1):i5(i6));
                Rz3 = Rz2(i5(1):i5(i6))';
            end
            i7 = find(diff(Rz3)>0);
            if(~isempty(i7))
                t_eye1 = t(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                Rz2 = Rz1(i3(1)-1:i3(1)-1+i5(1)+i7(1));
                i5 = find(Rz2<lpRr*min(Rz2) & Rz2>upRr*min(Rz2));
                i6 = find(diff(i5)>1);
                if(isempty(i6))
                    t_eye2 = t_eye1(i5(1):i5(end));
                    Rz3 = Rz2(i5(1):i5(end))';
                else
                    t_eye2 = t_eye1(i5(1):i5(i6));
                    Rz3 = Rz2(i5(1):i5(i6))';
                end
            end
            t_eye3 = [ones(length(t_eye2),1),t_eye2];
            b2 = t_eye3\Rz3;
            eyefit = t_chair3*b2;
            %Gain and latency calculations
            GaRr(i,1) = b2(2)/b(2); %ratio
            latRr(i,1) = 1000*(-b2(1)/b2(2) + b(1)/b(2)); %in ms
            %Calculate gain the AM way without latency shift
            i_s = find(t==t_eye2(1));
            AMGaRr(i,1) = mean(eyefit(i_s:i_s+length(Rz3)-1)./chairfit(i_s:i_s+length(Rz3)-1));
            %Plot
            if(plot_check)
                plot(t,chair1,'k',t,Rz1,'m',t,chairfit,'k',t,eyefit,'m',t(i_s:i_s+length(Rz3)-1),eyefit(i_s:i_s+length(Rz3)-1),'m*')
                title(['Check Fits for Right Eye Going to the Right on ',mouse])
                xlabel('Time (s)')
                ylabel('Angular Velocity (dps)')
                legend('Inverted Chair','Right Eye','Chair Fit','Right Eye Fit')
                axis([t(1) t(end) -1.3*maxvel 0.5*maxvel])
                pause;
            end
        end
        m_GaRr = mean(GaRr);
        m_latRr = mean(latRr);
        b3 = [m_GaRr*b(2)*(b(1)/b(2)-m_latRr/1000);m_GaRr*b(2)];
        m_eyefit = t_chair3*b3;        
        p1 = plot(t,chair1,'k');
        hold on
        p2 = plot(t,EyeRr,'m');
        p3 = plot(t,chairfit,'k','LineWidth',2);
        p4 = plot(t,m_eyefit,'m','LineWidth',2);
        title(['Check Fits for Right Eye Going to the Right on ',mouse])
        xlabel('Time (s)')
        ylabel('Angular Velocity (dps)')
        legend([p1,p2(1),p3,p4],{'Inverted Chair','Right Eye','Chair Fit','Right Eye Fit'})
        axis([t(1) t(end) -1.3*maxvel 0.5*maxvel])
        hold off       
        %Check if need to rerun at all
        rerun = input('Is the fit appropriate? (y/n): ','s');
        while ~(strcmp(rerun,'y')||strcmp(rerun,'n')) 
            disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
            rerun = input('Is the fit appropriate? (y/n): ','s');
        end
        if(strcmp(rerun,'n'))
            %Check the upper linear level
            u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            while ~(strcmp(u_rang,'y')||strcmp(u_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                u_rang = input('Is the limit on the linear range closest to the maximum velocity appropriate? (y/n): ','s');
            end
            if(strcmp(u_rang,'n'))
                u_str = ['The current upper limit percentage is ',num2str(upRr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upRr = new_rang/100;
            end
            %Check the lower linear level
            l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            while ~(strcmp(l_rang,'y')||strcmp(l_rang,'n')) 
                disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
                l_rang = input('Is the limit on the linear range closest to the stimulus onset appropriate? (y/n): ','s');
            end
            if(strcmp(l_rang,'n'))
                u_str = ['The current lower limit percentage is ',num2str(lpRr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpRr = new_rang/100;
            end
            %Make sure the limits are valid
            while(lpRr>=upRr)
                disp('The lower percentage limit cannot be greater than or equal to the upper percentage limit. They must be changed.')
                u_str = ['The current upper limit percentage is ',num2str(upRr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new upper limit percentage be? (range: 0-100) ');
                end
                upRr = new_rang/100;
                u_str = ['The current lower limit percentage is ',num2str(lpRr*100),' .'];
                disp(u_str)
                new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                while (new_rang<0||new_rang>100)
                    disp('Number out of range. Enter a pecentage between 0 and 100.')
                    new_rang = input('What should the new lower limit percentage be? (range: 0-100) ');
                end
                lpRr = new_rang/100;
            end
        else
            u_rang = 'y';
            l_rang = 'y';
        end
    end
    Rr.chair_m = chair1;
    Rr.chairfit = chairfit;
    Rr.eyefit = m_eyefit;
else
    disp('No traces accepted for right eye going right.')
    GaRr = [];
    latRr = [];
    AMGaRr = [];
    Rr.chair_m = NaN(1,length(t));
    Rr.chairfit = NaN(1,length(t));
    Rr.eyefit = NaN(1,length(t));
end
%% Make summary figures
mid = floor(length(t)/2);
subplot(2,1,1)
hold on
p1 = plot(t,Ll.chair_clean,'k');
plot(t,Rl.chair_clean,'k')
p3 = plot(t,Ll.chairfit,'k','LineWidth',2);
plot(t,Rl.chairfit,'k','LineWidth',2)
p5 = plot(t,Ll.eye_clean,'r');
p6 = plot(t,Ll.eyefit,'r','LineWidth',2);
p7 = plot(t,Rl.eye_clean,'m');
p8 = plot(t,Rl.eyefit,'m','LineWidth',2);
legend([p1(1),p3,p5(1),p6,p7(1),p8],{'Inverted Chair','Chair Fit','Left Eye','Left Eye Fit','Right Eye','Right Eye Fit'})
title('Leftward Eye Movements')
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
axis([t(1) t(mid) -0.5*maxvel 1.3*maxvel])
hold off
subplot(2,1,2)
hold on
p1 = plot(t,Lr.chair_clean,'k');
plot(t,Rr.chair_clean,'k')
p3 = plot(t,Lr.chairfit,'k','LineWidth',2);
plot(t,Rr.chairfit,'k','LineWidth',2)
p5 = plot(t,Lr.eye_clean,'r');
p6 = plot(t,Lr.eyefit,'r','LineWidth',2);
p7 = plot(t,Rr.eye_clean,'m');
p8 = plot(t,Rr.eyefit,'m','LineWidth',2);
legend([p1(1),p3,p5(1),p6,p7(1),p8],{'Inverted Chair','Chair Fit','Left Eye','Left Eye Fit','Right Eye','Right Eye Fit'})
title('Rightward Eye Movements')
xlabel('Time (s)')
ylabel('Angular Velocity (dps)')
axis([t(1) t(mid) -1.3*maxvel 0.5*maxvel])
suptitle([mouse,' Impulse Response'])
hold off
pause;
%% Make a summary table of everything from the Trace by Trace Analysis
meanGa = [mean(GaLl);mean(GaLr);mean(GaRl);mean(GaRr)];
stdGa = [std(GaLl);std(GaLr);std(GaRl);std(GaRr)];
meanlat = [mean(latLl);mean(latLr);mean(latRl);mean(latRr)];
stdlat = [std(latLl);std(latLr);std(latRl);std(latRr)];
meanAMGa = [mean(AMGaLl);mean(AMGaLr);mean(AMGaRl);mean(AMGaRr)];
stdAMGa = [std(AMGaLl);std(AMGaLr);std(AMGaRl);std(AMGaRr)];
n = [length(GaLl);length(GaLr);length(GaRl);length(GaRr)];
combo = zeros(1,7);
combo(1,1) = mean(meanGa,'omitnan');
combo(1,2) = sqrt(sum(stdGa.^2.*(n-1),'omitnan')/(sum(n,'omitnan')-length(stdGa(~isnan(stdGa)))));
combo(1,3) = mean(meanlat,'omitnan');
combo(1,4) = sqrt(sum(stdlat.^2.*(n-1),'omitnan')/(sum(n,'omitnan')-length(stdlat(~isnan(stdlat)))));
combo(1,5) = mean(meanAMGa,'omitnan');
combo(1,6) = sqrt(sum(stdAMGa.^2.*(n-1),'omitnan')/(sum(n,'omitnan')-length(stdAMGa(~isnan(stdAMGa)))));
combo(1,7) = sum(n,'omitnan');
summary_table = array2table([[meanGa,stdGa,meanlat,stdlat,meanAMGa,stdAMGa,n];combo]);
summary_table.Properties.RowNames = {'LeftEyeGoingLeft','LeftEyeGoingRight','RightEyeGoingLeft','RightEyeGoingRight','Combo'};
summary_table.Properties.VariableNames = {'meanGa','stdGa','meanlat','stdlat','meanAmericoGa','stdAmericoGa','n_cycles'};
params.mouse = mouse;
params.infile = info.fname;
params.analysis_time = clock;
params.minvel = minvel;
params.maxvel = maxvel;
params.fname = [mouse,'ImpulseDataAnalyzed.mat'];
params.lpLl = lpLl;
params.upLl = upLl;
params.lpLr = lpLr;
params.upLr = upLr;
params.lpRl = lpRl;
params.upRl = upRl;
params.lpRr = lpRr;
params.upRr = upRr;
Ll.GaLl = GaLl; %left eye going left
Lr.GaLr = GaLr; %left eye going right
Rl.GaRl = GaRl; %right eye going left
Rr.GaRr = GaRr; %right eye going right
Ll.latLl = latLl; %left eye going left
Lr.latLr = latLr; %left eye going right
Rl.latRl = latRl; %right eye going left
Rr.latRr = latRr; %right eye going right
Ll.AMGaLl = AMGaLl; %left eye going left
Lr.AMGaLr = AMGaLr; %left eye going right
Rl.AMGaRl = AMGaRl; %right eye going left
Rr.AMGaRr = AMGaRr; %right eye going right
ImpulseAnalysisOutput.t = t;
ImpulseAnalysisOutput.Ll = Ll;
ImpulseAnalysisOutput.Lr = Lr;
ImpulseAnalysisOutput.Rl = Rl;
ImpulseAnalysisOutput.Rr = Rr;
ImpulseAnalysisOutput.summary = summary_table;
ImpulseAnalysisOutput.info = params;
ImpulseAnalysisOutput.info.CleanImpulseData = info;
end