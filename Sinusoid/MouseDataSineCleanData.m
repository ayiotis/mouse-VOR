%% Mouse Data Sine Clean Data
%This function takes in the strcut SegData created from the Mouse Data Sine
%Segmentation code to remove erroneous traces and remove quick phases.
%It outputs a struct CleanDat with the erroneous cycles and quick phases
%removed.

%Written by Andrianna Ayiotis
%Last updated 06/22/2018

function CleanDat = MouseDataSineCleanData(SegDat)
%% Load data from file
t = SegDat.t;
cdats_seg = SegDat.Chair;
ldats_seg = SegDat.LEye;
rdats_seg = SegDat.REye;
freq = SegDat.info.true_freq;
round_freq = SegDat.info.round_freq;
mouse = SegDat.info.mouse;
maxvel = SegDat.info.maxvel;
%% Show plot of data to be cleaned first
plot_title = [num2str(mouse),' ',num2str(round_freq),' Hz Segmented Data'];
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
pause;
%% Delete any eroneous data traces 
%Make an average chair graph
m_chair = mean(cdats_seg,1);
fit = @(C_p,tt_seg)  C_p(1).*(sin(2*pi*C_p(2)*tt_seg) + C_p(3)*pi/180) + C_p(4);
%Definie the least squares cost function to minimize
LSCF = @(C_p) sum((fit(C_p,t) - m_chair).^2);                             
chair_params = fminsearch(LSCF, [(max(m_chair)-min(m_chair))/2;freq;0;0]);
chair_amp = chair_params(1);
chair_freq = chair_params(2);
Chair.All = cdats_seg;
Chair.cycle_avg = m_chair;
Chair.amp = chair_amp;
Chair.freq = chair_freq;
Chair.phase = chair_params(3);
Chair.offset = chair_params(4);
%Cycle through the traces and have the user choose to keep the trace or
%delete it for the left eye data
cycles = size(cdats_seg,1);
commit = false;
maxvel = Chair.amp*1.5;
a_thresh_l = 1000; %for quick phase acceleration threshold
qp_rm_l = false;
qp_rm_r = false;
while(~commit)
    rm_trace_l = [];
    chair_seg = cdats_seg;
    Leye_seg = ldats_seg;    
    for i = 1:cycles
        message = ['Trace ',num2str(i),' of ',num2str(cycles)];
        disp(message)
        p1 = plot(t,chair_seg,'k','LineWidth',5);
        hold on;
        p4 = plot(t,rdats_seg,'m');
        p2 = plot(t,Leye_seg,'r');
        p3 = plot(t,ldats_seg(i,:),'b','LineWidth',5);
        hold off;
        title('Check the trace quality')
        axis([0 t(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1),p2(1),p4(1),p3(1)],{'Chair','Left Eye','Right Eye','Trace'})
        u_keep = input('Keep this trace? (y/n)','s');
        while ~(strcmp(u_keep,'y')||strcmp(u_keep,'n')) 
           disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
           u_keep = input('Keep this trace? (y/n)','s');
        end
        if strcmp(u_keep,'n')
            Leye_seg(i-length(rm_trace_l),:) = [];
            chair_seg(i-length(rm_trace_l),:) = [];
            rm_trace_l = [rm_trace_l,i];
        end
    end    
    if(~isempty(Leye_seg))
        p1 = plot(t,chair_seg,'k');
        hold on;
        p3 = plot(t,rdats_seg,'m');
        p2 = plot(t,Leye_seg,'r');
        hold off;
        title('Check the new data quality')
        axis([0 t(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1) p2(1) p3(1)],{'Chair','Left Eye','Right Eye'})
        %Remove quick phases if needed
        rm_qp = input('Remove quick phases? (y/n) ','s');
        while(~strcmp(rm_qp,'y')&&~strcmp(rm_qp,'n'))
            disp('Invalid character entered. Only enter "y" or "n" without quotation marks.')
            rm_qp = input('Remove quick phases? (y/n) ','s');
        end
        while(strcmp(rm_qp,'y'))
            qp_rm_l = true;
            Lqpr = Leye_seg;
            Leye_a = diff(Leye_seg,1,2)/mean(diff(t));
            Leye_a_sized = [zeros(size(Leye_seg,1),1),Leye_a];
            Lqpr(abs(Leye_a_sized)>abs(a_thresh_l)) = nan;
            for j = 1:size(Lqpr,1)
                a_trace = Leye_a_sized(j,:)';
                nan_trace = Lqpr(j,:);
                inds = find(isnan(nan_trace))';
                if(~isempty(inds))
                    inds2 = [0;find(diff(inds)>1);length(inds)];
                    for k = 1:length(inds2)-1
                        i1 = inds(inds2(k)+1);
                        i2 = inds(inds2(k+1));
                        front_snip = a_trace(i1-1:-1:1);
                        back_snip = a_trace(i2+1:end);
                        %Test whether the nan area is increasing or decreasing 
                        sign = mean(a_trace(i1:i2))>0;
                        if(sign)
                            i3 = length(front_snip)-find(front_snip<=0,1)+1;
                            i4 = i2 + find(back_snip<=0,1);
                        else
                            i3 = length(front_snip)-find(front_snip>=0,1)+1;
                            i4 = i2 + find(back_snip>=0,1);
                        end
                        Lqpr(j,i3:i4) = nan;
                    end
                end
            end
            %Plot the results of QP removal
            p1 = plot(t,chair_seg,'k');
            hold on;
            p4 = plot(t,rdats_seg,'m');
            p2 = plot(t,Leye_seg,'r');
            p3 = plot(t,Lqpr,'b');
            hold off
            legend([p1(1),p2(1),p4(1),p3(1)],{'Chair','Left Eye','Right Eye','QSP gone'})
            axis([0 t(end) -1.1*maxvel 1.1*maxvel])
            %Try again if needed
            message = ['The current acceleration threshold for a saccade is ',num2str(a_thresh_l), 'deg/s^2.'];
            disp(message)
            rm_qp = input('Would you like to change this value? (y/n)','s');
            while(~strcmp(rm_qp,'y')&&~strcmp(rm_qp,'n'))
                disp('Invalid character entered. Only enter "y" or "n" without quotation marks.')
                rm_qp = input('Would you like to change this value? (y/n)','s');
            end
            if(strcmp(rm_qp,'y'))
                a_thresh_l = input('Enter new value: ');
            else
                Leye_seg = Lqpr;
            end
        end
    else
        disp('No traces accepted')
        Leye_seg = NaN(1,length(t));
    end    
    message = [num2str(cycles-length(rm_trace_l)),' of ',num2str(cycles),' total traces accepted.'];
    disp(message)
    u_commit = input('Commit changes (y) or try the whole process again (n)?','s');
    while ~(strcmp(u_commit,'y')||strcmp(u_commit,'n')) 
       disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
       u_commit = input('Commit changes (y) or try the whole process again (n)?','s');
    end
    if strcmp(u_commit,'y')
       commit = true;
    end    
end

%Cycle through the traces and have the user choose to keep the trace or
%remove it for the right eye data
a_thresh_r = a_thresh_l;
commit = false;
while(~commit)
    rm_trace_r = [];
    chair_seg = cdats_seg;
    Reye_seg = rdats_seg;
    for i = 1:cycles
        message = ['Trace ',num2str(i),' of ',num2str(cycles)];
        disp(message)
        p1 = plot(t,chair_seg,'k','LineWidth',5);
        hold on;
        p4 = plot(t,Leye_seg,'r');
        p2 = plot(t,Reye_seg,'m');
        p3 = plot(t,rdats_seg(i,:),'b','LineWidth',5);
        hold off;
        title('Check the trace quality')
        axis([0 t(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1) p4(1) p2(1) p3(1)],{'Chair','Left Eye','Right Eye','Trace'})
        u_keep = input('Keep this trace? (y/n)','s');
        while ~(strcmp(u_keep,'y')||strcmp(u_keep,'n')) 
           disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
           u_keep = input('Keep this trace? (y/n)','s');
        end
        if strcmp(u_keep,'n')
            chair_seg(i-length(rm_trace_r),:) = [];
            Reye_seg(i-length(rm_trace_r),:) = [];
            rm_trace_r = [rm_trace_r,i];
        end
    end
    
    if(~isempty(Reye_seg))
        p1 = plot(t,chair_seg,'k');
        hold on;
        p3 = plot(t,Leye_seg,'r');
        p2 = plot(t,Reye_seg,'m');
        hold off;
        title('Check the new data quality')
        axis([0 t(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1) p2(1) p3(1)],{'Chair','Left Eye','Right Eye'})        
        %Remove quick phases if needed
        rm_qp = input('Remove quick phases? (y/n) ','s');
        while(~strcmp(rm_qp,'y')&&~strcmp(rm_qp,'n'))
            disp('Invalid character entered. Only enter "y" or "n" without quotation marks.')
            rm_qp = input('Remove quick phases? (y/n) ','s');
        end
        while(strcmp(rm_qp,'y'))  
            qp_rm_r = true;
            Rqpr = Reye_seg;
            Reye_a = diff(Reye_seg,1,2)/mean(diff(t));
            Reye_a_sized = [zeros(size(Reye_seg,1),1),Reye_a];
            Rqpr(abs(Reye_a_sized)>abs(a_thresh_r)) = nan;
            for j = 1:size(Rqpr,1)
                a_trace = Reye_a_sized(j,:)';
                nan_trace = Rqpr(j,:);
                inds = find(isnan(nan_trace))';
                if(~isempty(inds))
                    inds2 = [0;find(diff(inds)>1);length(inds)];
                    for k = 1:length(inds2)-1
                        i1 = inds(inds2(k)+1);
                        i2 = inds(inds2(k+1));
                        front_snip = a_trace(i1-1:-1:1);
                        back_snip = a_trace(i2+1:end);
                        %Test whether the nan area is increasing or decreasing 
                        sign = mean(a_trace(i1:i2))>0;
                        if(sign)
                            i3 = length(front_snip)-find(front_snip<=0,1)+1;
                            i4 = i2 + find(back_snip<=0,1);
                        else
                            i3 = length(front_snip)-find(front_snip>=0,1)+1;
                            i4 = i2 + find(back_snip>=0,1);
                        end
                        Rqpr(j,i3:i4) = nan;
                    end
                end
            end
            %Plot the results of QP removal
            p1 = plot(t,chair_seg,'k');
            hold on;
            p4 = plot(t,Leye_seg,'r');
            p2 = plot(t,Reye_seg,'m');
            p3 = plot(t,Rqpr,'b');
            hold off
            legend([p1(1),p2(1),p4(1),p3(1)],{'Chair','Left Eye','Right Eye','QSP gone'})
            axis([0 t(end) -1.1*maxvel 1.1*maxvel])            
            %Try again if needed
            message = ['The current acceleration threshold for a saccade is ',num2str(a_thresh_r), 'deg/s^2.'];
            disp(message)
            rm_qp = input('Would you like to change this value? (y/n)','s');
            while(~strcmp(rm_qp,'y')&&~strcmp(rm_qp,'n'))
                disp('Invalid character entered. Only enter "y" or "n" without quotation marks.')
                rm_qp = input('Would you like to change this value? (y/n)','s');
            end
            if(strcmp(rm_qp,'y'))
                a_thresh_r = input('Enter new value: ');
            else
                Reye_seg = Rqpr;
            end
        end
    else
        disp('No traces accepted')
        Reye_seg = NaN(1,length(t));
    end    
    message = [num2str(cycles-length(rm_trace_r)),' of ',num2str(cycles),' total traces accepted.'];
    disp(message)
    u_commit = input('Commit changes (y) or try the process again (n)?','s');
    while ~(strcmp(u_commit,'y')||strcmp(u_commit,'n')) 
       disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
       u_commit = input('Commit changes (y) or try the process again (n)?','s');
    end
    if strcmp(u_commit,'y')
       commit = true; 
    end    
end
%% Plot the accepted traces
p1 = plot(t,m_chair,'k','LineWidth',2);
hold on
p2 = plot(t,Leye_seg,'Color','r'); 
p3 = plot(t,Reye_seg,'Color','m');
hold off
legend([p1,p2(1),p3(1)],{'Chair','Left Eye Traces','Right Eye Traces'})
long_title = [mouse,' ',num2str(round_freq),' Hz Cycles'];
title(long_title)
xlabel('Time (s)')
ylabel('Velocity (dps)')
pause;
%% Compile data
info = SegDat.info;
info.SegDat = SegDat.info;
info = rmfield(info,'seg_num');
info.infile = info.fname;
info.analysis_time = clock;
info.rm_qp = qp_rm_l|qp_rm_r;
info.qp_a_thresh_l = a_thresh_l;
info.qp_a_thresh_r = a_thresh_r;
fname = [info.fname(1:end-4),'Clean.mat'];
info.fname = fname;
L.All = ldats_seg;
L.bad_traces = rm_trace_l;
L.cycles = Leye_seg;
R.All = rdats_seg;
R.bad_traces = rm_trace_r;
R.cycles = Reye_seg;
CleanDat.info = info;
CleanDat.t = t;
CleanDat.Chair = Chair;
CleanDat.LEye = L;
CleanDat.REye = R;
end