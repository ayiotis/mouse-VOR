%% Mouse Data Sine Segmentation
%This function takes in the mat files created from EMA, segments them into
%different frequencies and then alligns them for further processing.
%It takes in the arguments of mouse ID, a struct called SineData, a boolean
%that allows the user to save each of the segments of the file, and a 
%boolean that decides whether to preview the data during segmentation.
%It outputs a struct SegDat that contains the segmented data.
%The input file must have a struct called SineData with fields Time, Chair,
%Lz, and Rz.

%Written by Andrianna Ayiotis
%Last updated 06/13/2018

function MouseDataSineSegmentation(mouse,in_fname,SineData,save_data,prev_data,check_align)
disp(mouse)
disp(in_fname)
%% Get data from file
time = SineData.Time;
time = time - time(1); 
chair = SineData.Chair;
LEye = SineData.Lz;
REye = SineData.Rz;
maxvel = max(chair);
%% Segment into different frequencies using the chair velocity
l = length(time);
% This accounts for whether there is zero padding at the beginning of the
% data or not
if chair(1) == 0
    starts = 1;
else
   starts = []; 
end
stops = [];
%Find the string of 0s in the array
for i = 2:l
    if (chair(i) == 0 && chair(i-1) ~= 0) %Beginning of a sequence of 0s
        starts = [starts;i];
    elseif (chair(i) ~= 0 && chair(i-1) == 0) %End of a sequence of 0s
        stops = [stops;i-1];
    end
end
%%Accounts for 0 padding at the end of an array if needed
if length(starts) == length(stops)
    res = [starts,stops,stops-starts+1];
elseif length(starts) == length(stops)+1
    stops = [stops;l];
    res = [starts,stops,stops-starts+1];
else
    disp('Discrepancy between number of start and stop points.')
end
%Find only the places that have more than 100 consecutive 0s
breaks = res((res(:,3)>100),:);
blen = size(breaks,1);
%From the breaks, I can find the segments of data
%If no 0 padding at start
if breaks(1,1)-100 > 0
    segs(1,1) = 1;
    segs(1,2) = breaks(1,1)-1;
else
    segs = [];
end
for i = 1:blen-1
    segs = [segs; breaks(i,2)+1,breaks(i+1,1)-1];
end
%If no 0 padding at end
if l - breaks(blen,2) > 100
    segs = [segs; breaks(blen,2)+1,l];
end
%Make sure the cycle starts on a 0
segs(:,1) = segs(:,1)-1;
if(segs(1,1)==0)
   segs(1,1) = 1; 
end
p1 = plot(time,chair,'k');
hold on
p2 = plot(time,LEye,'r');
p3 = plot(time,REye,'m');
hold off
for i = 1:size(segs,1)
    h1 = line([time(segs(i,1)) time(segs(i,1))],[-1.1*maxvel 1.1*maxvel]);
    h2 = line([time(segs(i,2)) time(segs(i,2))],[-1.1*maxvel 1.1*maxvel]);
    % Set properties of lines
    set([h1 h2],'Color','k','LineWidth',2)
    % Add a patch
    patch([time(segs(i,1)) time(segs(i,2)) time(segs(i,2)) time(segs(i,1))],[-1.1*maxvel -1.1*maxvel 1.1*maxvel 1.1*maxvel],'g')
end
set(gca,'children',flipud(get(gca,'children')))
title('Check segments of data in the file')
legend([p1,p2,p3],{'Chair','Left Eye','Right Eye'})
xlabel('Time (s)')
ylabel('Deg/s')
axis([0 time(end) -1.1*maxvel 1.1*maxvel])
pause;
%% Analyze the data separately and save each alligned segment
slen = size(segs,1);
%Ensure t vector is rounded correctly (artifact of data input)
dt = mean(diff(time));
time = (0:dt:(length(time)-1)*dt)';
for i = 1:slen
    t = time(segs(i,1):segs(i,2),1);
    cdat = chair(segs(i,1):segs(i,2),1);
    ldat = LEye(segs(i,1):segs(i,2),1);
    rdat = REye(segs(i,1):segs(i,2),1);
    % Interpolate data points to allign
    tt_c = t(1):0.001:t(end);
    cdats_c = spline(t,cdat,tt_c);
    ldats_c = spline(t,ldat,tt_c);
    rdats_c = spline(t,rdat,tt_c);
    %% Find the native frequency of the chair velocity
    Fs=1/0.001;
    N=length(cdats_c);
    fx=fft(cdats_c)/N;
    f=(0:N-1)*Fs/N;
    x = f(1:floor(end/2));
    y = 10*log10(abs(fx(1:floor(end/2))));
    [~,ind] = max(y);
    freq = x(ind);
    if(freq > 0.9) %This should get frequencies 10, 5, 2, and 1 Hz.
        round_freq = round(freq);
    elseif (freq>0.08&&freq<=0.9) %This should get 0.5, 0.2, and 0.1 Hz.
        round_freq = round(freq,1);
    else %This should get 0.05 and 0.02 Hz.
        round_freq = round(freq,2);
    end
    %% Trim and align the chair velocity traces
    cdats = cdats_c;
    ldats = ldats_c;
    rdats = rdats_c;
    tt = tt_c;    
    if(round_freq>=0.5)
        cdats = cdats(50:end-50);
        ldats = ldats(50:end-50);
        rdats = rdats(50:end-50);
        tt = tt(50:end);
    end    
    if cdats(1) < 0
       ind1 = find(cdats>=0,1);
    else
       ind0 = find(cdats<=0,1);
       tt = tt(ind0+1:end);
       cdats = cdats(ind0+1:end);
       ldats = ldats(ind0+1:end);
       rdats = rdats(ind0+1:end);
       ind1 = find(cdats>=0,1);
    end    
    %Using the rounded frequency is normally right
    points = round(1/(round_freq*0.001));
    %Trim the beginning to start at 0
    tt = tt(ind1:end);
    cdats = cdats(ind1:end);
    ldats = ldats(ind1:end);
    rdats = rdats(ind1:end);
    tt = tt-tt(1);
    %Trim the end to be a good stopping point
    cycles = floor(length(cdats)/points);
    if(round_freq>=5)
        cycles = cycles - 1;
    end    
    tt = tt(1:cycles*points);
    cdats = cdats(1:cycles*points);
    ldats = ldats(1:cycles*points);
    rdats = rdats(1:cycles*points);
    tt_seg = tt(1:points);
    cdats_seg = reshape(cdats,points,cycles)';
    %Confirm allignment worked visually
    plot(tt_seg,cdats_seg,'k')
    title('Confirm proper alignment of stimulus traces')
    xlabel('Time (s)')
    ylabel('Chair vel (dps)')
    axis([0 tt_seg(end) -1.1*maxvel 1.1*maxvel])
    %Check that there are the right number of points
    if(check_align)
        g_align = input('Are the traces properly aligned? (y/n) ','s');
        while(~strcmp(g_align,'n')&&~strcmp(g_align,'y'))
            disp('Invalid character entered. Please enter only "y" or "n" without quotation marks.')
            g_align = input('Are the traces properly aligned?(y/n) ','s');
        end
    else
        g_align = 'y';
    end
    while(strcmp(g_align,'n'))
        message = ['There are currently ',num2str(points),' points per cycle.'];
        disp(message)
        points = input('How many points would be appropriate? ');
        cdats = cdats_c;
        ldats = ldats_c;
        rdats = rdats_c;
        tt = tt_c;       
        if(round_freq>=0.5)
            cdats = cdats(50:end-50);
            ldats = ldats(50:end-50);
            rdats = rdats(50:end-50);
            tt = tt(50:end);
        end       
        if cdats(1) < 0
           ind1 = find(cdats>=0,1);
        else
           ind0 = find(cdats<=0,1);
           tt = tt(ind0+1:end);
           cdats = cdats(ind0+1:end);
           ldats = ldats(ind0+1:end);
           rdats = rdats(ind0+1:end);
           ind1 = find(cdats>=0,1);
        end
        %Trim the beginning to start at 0
        tt = tt(ind1:end);
        cdats = cdats(ind1:end);
        ldats = ldats(ind1:end);
        rdats = rdats(ind1:end);
        tt = tt-tt(1);    
        %Trim the end to be a good stopping point
        cycles = floor(length(cdats)/points);
        if(round_freq>=5)
            cycles = cycles - 1;
        end        
        tt = tt(1:cycles*points);
        cdats = cdats(1:cycles*points);
        ldats = ldats(1:cycles*points);
        rdats = rdats(1:cycles*points);
        tt_seg = tt(1:points);
        cdats_seg = reshape(cdats,points,cycles)';
        %Confirm alignment worked visually
        plot(tt_seg,cdats_seg,'k')
        title('Confirm proper alignment of stimulus traces')
        xlabel('Time (s)')
        ylabel('Chair vel (dps)')
        axis([0 tt_seg(end) -1.1*maxvel 1.1*maxvel])
        g_align = input('Are the traces properly aligned? (y/n) ','s');
        while(~strcmp(g_align,'n')&&~strcmp(g_align,'y'))
            disp('Invalid character entered. Please enter only "y" or "n" without quotation marks.')
            g_align = input('Are the traces properly aligned?(y/n) ','s');
        end   
    end
    ldats_seg = reshape(ldats,points,cycles)';
    rdats_seg = reshape(rdats,points,cycles)';    
    if(prev_data)
        plot_title = [num2str(round_freq),' Hz Segmented Data'];
        %Preview the data
        p1 = plot(tt_seg,cdats_seg,'k');
        hold on
        p2 = plot(tt_seg,ldats_seg,'r');
        p3 = plot(tt_seg,rdats_seg,'m');
        hold off
        title(plot_title)
        xlabel('Time (s)')
        ylabel('Velocity (dps)')
        axis([0 tt_seg(end) -1.1*maxvel 1.1*maxvel]) 
        legend([p1(1),p2(1),p3(1)],{'Chair','Left Eye','Right Eye'})
        pause;
    end
    %Catch multiple versions (up to 3)
    out_file = [mouse,'-',strrep(num2str(round_freq),'.','p'),'HzSineSegmented.mat'];
    if(exist(out_file,'file'))
        out_file = [mouse,'-',strrep(num2str(round_freq),'.','p'),'HzSineSegmented_v2.mat'];
        if(exist(out_file,'file'))
            out_file = [mouse,'-',strrep(num2str(round_freq),'.','p'),'HzSineSegmented_v3.mat'];
            disp('Now three versions of this file exist. More versions will write over v3.')
        end
    end    
    %Make struct to save
    info.infile = in_fname;
    info.mouse = mouse;
    info.analysis_time = clock;
    info.round_freq = round_freq;
    info.true_freq = freq;
    info.seg_num = i;
    info.maxvel = maxvel;
    info.points_p_cyle = points;
    info.fname = out_file;

    SegDat.t = tt_seg;
    SegDat.Chair = cdats_seg;
    SegDat.LEye = ldats_seg;
    SegDat.REye = rdats_seg;
    SegDat.info = info;   
    
    if(save_data)
       save(out_file,'SegDat') 
    end
end
end