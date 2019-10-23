%% Mouse Data Impulse Clean Data
% This function takes in the mouseID string, the name of the input file and
% a table called ImpulseData with columns "Time", "Lz", "Rz", and "Chair"
% (extra columns are allowed but will not be analyzed).
% The data are segmented and the user selects which traces (if any) should
% be discarded. 
% The function outputs a struct called CleanImpulseData

%Made by Andrianna Ayiotis
%Last updated 08/08/2018
function CleanImpulseData = MouseDataImpulseCleanData(mouse,infile,ImpulseData)
%% Load needed variables
time = ImpulseData.Time;
Lz = ImpulseData.Lz;
Rz = ImpulseData.Rz;
Chair = ImpulseData.Chair;
%User defined parameters--change manually if needed
minvel = 60; %dps (minimum chair velocity for alignment)
maxvel = 300; %dps (maximum chair velocity for alignment)
dur = 700; %ms (duration of plot after the minvel is reached)
pre_imp = 150; %ms (how many ms on the plot before the minvel is reached)
%% Segment the Chair Velocity
%Trim the ends 
i1 = find(Chair == 0);
off = i1(1);
Chair_t = Chair(off:i1(end));
%Movements to the left
i2 = find(Chair_t > minvel & Chair_t < 1.1*maxvel);
i3 = find(diff(i2)>1);
starts = [i2(1);i2(i3+1)];
ms_p_ind = 1000*(time(end) - time(1))/length(time); %milliseconds per index
indf = floor(pre_imp/ms_p_ind);
inda = ceil(dur/ms_p_ind);
starts_off = starts - 1 - indf + off;
ends_off = starts - 1 + off + inda;
for i = 1:length(starts_off)
    chair_seg_l(i,:) = Chair(starts_off(i):ends_off(i))';
    Lz_seg_l(i,:) = Lz(starts_off(i):ends_off(i))';
    Rz_seg_l(i,:) = Rz(starts_off(i):ends_off(i))';
end
t_seg_l = time(1:size(chair_seg_l,2));
%Movements to the right 
i2 = find(Chair_t < -minvel & Chair_t > -1.1*maxvel);
i3 = find(diff(i2)>1);
starts = [i2(1);i2(i3+1)];
ms_p_ind = 1000*(time(end) - time(1))/length(time); %milliseconds per index
indf = floor(pre_imp/ms_p_ind);
inda = ceil(dur/ms_p_ind);
starts_off = starts - 1 - indf + off;
ends_off = starts - 1 + off + inda;
for i = 1:length(starts_off)
    chair_seg_r(i,:) = Chair(starts_off(i):ends_off(i))';
    Lz_seg_r(i,:) = Lz(starts_off(i):ends_off(i))';
    Rz_seg_r(i,:) = Rz(starts_off(i):ends_off(i))';
end
t_seg_r = time(1:size(chair_seg_r,2));
subplot(2,1,1)
plot(t_seg_l,chair_seg_l,'k')
hold on
plot(t_seg_l,Lz_seg_l,'r')
plot(t_seg_l,Rz_seg_l,'m')
title('Leftward Eye Movements')
ylabel('Angular Velocity (dps)')
axis([t_seg_l(1) t_seg_l(end) -1.1*maxvel 1.1*maxvel])
hold off
subplot(2,1,2)
p1 = plot(t_seg_r,chair_seg_r,'k');
hold on
p2 = plot(t_seg_r,Lz_seg_r,'r');
p3 = plot(t_seg_r,Rz_seg_r,'m');
title('Rightward Eye Movements')
ylabel('Angular Velocity (dps)')
xlabel('Time (s)')
axis([t_seg_r(1) t_seg_r(end) -1.1*maxvel 1.1*maxvel])
legend([p1(1) p2(1) p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
suptitle(['Check segmenting for ',mouse])
hold off
pause;
close;
%% Remove bad traces (Left eye going to the left)
commit = 'n';
while(strcmp(commit,'n'))
    rm_trace_Ll = [];
    chairsegl_Lc = chair_seg_l;
    Lzsegl_c = Lz_seg_l;
    for i = 1:size(chair_seg_l,1)
        disp(['Trace ',num2str(i),' of ',num2str(size(chair_seg_l,1))])
        p1 = plot(t_seg_l,chair_seg_l,'k');
        hold on;
        p2 = plot(t_seg_l,Rz_seg_l,'m');
        p3 = plot(t_seg_l,Lzsegl_c,'r');
        p4 = plot(t_seg_l,Lz_seg_l(i,:),'b','LineWidth',5);
        title(['Check Left Eye Moving to the Left for ',mouse])
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        axis([t_seg_l(1) t_seg_l(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1),p3(1),p2(1),p4(1)],{'Inverted Chair','Left Eye','Right Eye','Trace'})
        hold off;
        u_keep = input('Keep this trace? (y/n)','s');
        while ~(strcmp(u_keep,'y')||strcmp(u_keep,'n')) 
            disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
            u_keep = input('Keep this trace? (y/n)','s');
        end
        if strcmp(u_keep,'n')
            chairsegl_Lc(i-length(rm_trace_Ll),:) = [];
            Lzsegl_c(i-length(rm_trace_Ll),:) = [];
            rm_trace_Ll = [rm_trace_Ll,i];
        end
    end
    if(~isempty(Lzsegl_c))
        p1 = plot(t_seg_l,chairsegl_Lc,'k');
        hold on;
        p3 = plot(t_seg_l,Rz_seg_l,'m');
        p2 = plot(t_seg_l,Lzsegl_c,'r');
        title(['Check Left Eye Moving to the Left for ',mouse])
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        axis([t_seg_l(1) t_seg_l(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1),p2(1),p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
        hold off;
    else
        disp('No traces accepted')
    end
    %Get user input on whether to repeat the process
    message = [num2str(size(chair_seg_l,1)-length(rm_trace_Ll)),' of ',num2str(size(chair_seg_l,1)),' total traces accepted.'];
    disp(message)
    commit = input('Commit changes (y) or try the process again (n)?','s');
    while ~(strcmp(commit,'y')||strcmp(commit,'n')) 
       disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
       commit = input('Commit changes (y) or try the process again (n)?','s');
    end
end
Ll.chair_all = chair_seg_l;
Ll.eye_all = Lz_seg_l;
if(~isempty(Lzsegl_c))
    Ll.chair_clean = chairsegl_Lc;
    Ll.eye_clean = Lzsegl_c; 
else
    Ll.chair_clean = NaN(1,length(t_seg_l));
    Ll.eye_clean = NaN(1,length(t_seg_l));
end
Ll.rm_trace = rm_trace_Ll;
%% Remove bad traces (Right eye going to the left)
commit = 'n';
while(strcmp(commit,'n'))
    rm_trace_Rl = [];
    chairsegl_Rc = chair_seg_l;
    Rzsegl_c = Rz_seg_l;
    for i = 1:size(chair_seg_l,1)
        disp(['Trace ',num2str(i),' of ',num2str(size(chair_seg_l,1))])
        p1 = plot(t_seg_l,chair_seg_l,'k');
        hold on;
        p2 = plot(t_seg_l,Lzsegl_c,'r');
        p3 = plot(t_seg_l,Rzsegl_c,'m');
        p4 = plot(t_seg_l,Rz_seg_l(i,:),'b','LineWidth',5);
        legend([p1(1),p2(1),p3(1),p4(1)],{'Inverted Chair','Left Eye','Right Eye','Trace'})   
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        title(['Check Right Eye Going to the Left for ',mouse])
        axis([t_seg_l(1) t_seg_l(end) -1.1*maxvel 1.1*maxvel])
        hold off;
        u_keep = input('Keep this trace? (y/n)','s');
        while ~(strcmp(u_keep,'y')||strcmp(u_keep,'n')) 
           disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
           u_keep = input('Keep this trace? (y/n)','s');
        end
        if strcmp(u_keep,'n')
            chairsegl_Rc(i-length(rm_trace_Rl),:) = [];
            Rzsegl_c(i-length(rm_trace_Rl),:) = [];
            rm_trace_Rl = [rm_trace_Rl,i];
        end
    end
    if(~isempty(Rzsegl_c))
        p1 = plot(t_seg_l,chairsegl_Rc,'k');
        hold on;
        p2 = plot(t_seg_l,Lzsegl_c,'r');
        p3 = plot(t_seg_l,Rzsegl_c,'m');
        title(['Check Right Eye Moving to the Left for ',mouse])
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        axis([t_seg_l(1) t_seg_l(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1) p2(1) p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
        hold off;
    end
    message = [num2str(size(chair_seg_l,1)-length(rm_trace_Rl)),' of ',num2str(size(chair_seg_l,1)),' total traces accepted.'];
    disp(message)
    commit = input('Commit changes (y) or try the process again (n)?','s');
    while ~(strcmp(commit,'y')||strcmp(commit,'n')) 
       disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
       commit = input('Commit changes (y) or try the process again (n)?','s');
    end
end
Rl.chair_all = chair_seg_l;
Rl.eye_all = Rz_seg_l;
if(~isempty(Rzsegl_c))
    Rl.chair_clean = chairsegl_Rc;
    Rl.eye_clean = Rzsegl_c; 
else
    Rl.chair_clean = NaN(1,length(t_seg_l));
    Rl.eye_clean = NaN(1,length(t_seg_l));
end
Rl.rm_trace = rm_trace_Rl;
%% Remove bad traces (Left eye going to the right)
commit = 'n';
while(strcmp(commit,'n'))
    rm_trace_Lr = [];
    chairsegr_Lc = chair_seg_r;
    Lzsegr_c = Lz_seg_r;
    for i = 1:size(chair_seg_r,1)
        message = ['Trace ',num2str(i),' of ',num2str(size(chair_seg_r,1))];
        disp(message)
        p1 = plot(t_seg_r,chair_seg_r,'k');
        hold on;
        p2 = plot(t_seg_r,Rz_seg_r,'m');
        p3 = plot(t_seg_r,Lzsegr_c,'r');
        p4 = plot(t_seg_r,Lz_seg_r(i,:),'b','LineWidth',5);
        title(['Check Left Eye Moving to the Right for ',mouse])
        axis([t_seg_r(1) t_seg_r(end) -1.1*maxvel 1.1*maxvel])
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        legend([p1(1),p3(1),p2(1),p4(1)],{'Inverted Chair','Left Eye','Right Eye','Trace'})
        hold off;
        u_keep = input('Keep this trace? (y/n)','s');
        while ~(strcmp(u_keep,'y')||strcmp(u_keep,'n')) 
           disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
           u_keep = input('Keep this trace? (y/n)','s');
        end
        if strcmp(u_keep,'n')
            chairsegr_Lc(i-length(rm_trace_Lr),:) = [];
            Lzsegr_c(i-length(rm_trace_Lr),:) = [];
            rm_trace_Lr = [rm_trace_Lr,i];
        end
    end
    if(~isempty(Lzsegr_c))
        p1 = plot(t_seg_r,chairsegr_Lc,'k');
        hold on;
        p3 = plot(t_seg_r,Rz_seg_r,'m');
        p2 = plot(t_seg_r,Lzsegr_c,'r');
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        title(['Check Left Eye Moving to the Right for ',mouse])
        axis([t_seg_r(1) t_seg_r(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1),p2(1),p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
        hold off;
    end
    message = [num2str(size(chair_seg_r,1)-length(rm_trace_Lr)),' of ',num2str(size(chair_seg_r,1)),' total traces accepted.'];
    disp(message)
    commit = input('Commit changes (y) or try the process again (n)?','s');
    while ~(strcmp(commit,'y')||strcmp(commit,'n')) 
       disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
       commit = input('Commit changes (y) or try the process again (n)?','s');
    end
end
Lr.chair_all = chair_seg_r;
Lr.eye_all = Lz_seg_r;
if(~isempty(Lzsegr_c))
    Lr.chair_clean = chairsegr_Lc;
    Lr.eye_clean = Lzsegr_c; 
else
    Lr.chair_clean = NaN(1,length(t_seg_r));
    Lr.eye_clean = NaN(1,length(t_seg_r));
end
Lr.rm_trace = rm_trace_Lr;
%% Remove bad traces (Right eye going to the right)
commit = 'n';
while(strcmp(commit,'n'))
    rm_trace_Rr = [];
    chairsegr_Rc = chair_seg_r;
    Rzsegr_c = Rz_seg_r;
    for i = 1:size(chair_seg_r,1)
        message = ['Trace ',num2str(i),' of ',num2str(size(chair_seg_r,1))];
        disp(message)
        p1 = plot(t_seg_r,chair_seg_r,'k');
        hold on;
        p2 = plot(t_seg_r,Lzsegr_c,'r');
        p3 = plot(t_seg_r,Rzsegr_c,'m');
        p4 = plot(t_seg_r,Rz_seg_r(i,:),'b','LineWidth',5);
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        legend([p1(1),p2(1),p3(1),p4(1)],{'Inverted Chair','Left Eye','Right Eye','Trace'})
        title(['Check Right Eye Moving to the Right for ',mouse])
        axis([t_seg_r(1) t_seg_r(end) -1.1*maxvel 1.1*maxvel])
        hold off;
        u_keep = input('Keep this trace? (y/n)','s');
        while ~(strcmp(u_keep,'y')||strcmp(u_keep,'n')) 
           disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
           u_keep = input('Keep this trace? (y/n)','s');
        end
        if strcmp(u_keep,'n')
            chairsegr_Rc(i-length(rm_trace_Rr),:) = [];
            Rzsegr_c(i-length(rm_trace_Rr),:) = [];
            rm_trace_Rr = [rm_trace_Rr,i];
        end
    end
    if(~isempty(Rzsegr_c))
        p1 = plot(t_seg_r,chairsegr_Rc,'k');
        hold on;
        p2 = plot(t_seg_r,Lzsegr_c,'r');
        p3 = plot(t_seg_r,Rzsegr_c,'m');
        title(['Check Right Eye Moving to the Right for ',mouse])
        ylabel('Angular Velocity (dps)')
        xlabel('Time (s)')
        axis([t_seg_r(1) t_seg_r(end) -1.1*maxvel 1.1*maxvel])
        legend([p1(1),p2(1),p3(1)],{'Inverted Chair','Left Eye','Right Eye'})
        hold off;
    end
    message = [num2str(size(chair_seg_r,1)-length(rm_trace_Rr)),' of ',num2str(size(chair_seg_r,1)),' total traces accepted.'];
    disp(message)
    commit = input('Commit changes (y) or try the process again (n)?','s');
    while ~(strcmp(commit,'y')||strcmp(commit,'n')) 
       disp('Invalid character entered. Only enter "y" for yes or "n" for no (without quotation marks).')
       commit = input('Commit changes (y) or try the process again (n)?','s');
    end
end
Rr.chair_all = chair_seg_r;
Rr.eye_all = Rz_seg_r;
if(~isempty(Rzsegr_c))
    Rr.chair_clean = chairsegr_Rc;
    Rr.eye_clean = Rzsegr_c; 
else
    Rr.chair_clean = NaN(1,length(t_seg_r));
    Rr.eye_clean = NaN(1,length(t_seg_r));
end
Rr.rm_trace = rm_trace_Rr;
%% Save it all in one struct for ease of use
info.mouse = mouse;
info.analysis_time = clock;
info.infile = infile;
info.fname = [mouse,'CleanImpulseData.mat'];
info.minvel = minvel;
info.maxvel = maxvel;
info.dur = dur;
info.pre_imp = pre_imp;
info.ms_p_ind = ms_p_ind;
CleanImpulseData.t = t_seg_l; %same as t_seg_r
CleanImpulseData.Ll = Ll;
CleanImpulseData.Rl = Rl;
CleanImpulseData.Lr = Lr;
CleanImpulseData.Rr = Rr;
CleanImpulseData.info = info;
end