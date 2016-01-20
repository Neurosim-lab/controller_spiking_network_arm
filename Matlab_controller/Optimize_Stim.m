%testscript
clear all
close all
clc

load('iterparams.mat');

%% params
load([dataFolder dataSubfolder 'Inv_Plant_target-',num2str(target_num),'_ptype-',num2str(perturbation_type),'_pperc-',num2str(perturbation_percentage)   '_' probing '.mat']);

M_Cell_Size = 192; 
cell_id = 101 - 1; %1:192
maxRate = 250; %max spiking frequency
stim_starting_time = 100; %stim starting time 100ms (or 400ms)

%% Data Formatting for Output (EM-Muscle) Response
%========================================================================================================================

whichfile = [bmmFolderOriginal,'/target-' num2str(target_num) '_ptype-1_pperc-0_cell-0_start-0_dur-0_rate-0_spk.mat'];
spk = load(whichfile);

%Partition into input segments (spike trains)
EM_muscles_spike_timing = {};
spike_count = length(spk.tspike); %total number of spikes in the set,

%Allocation spike timing array for each of the M_Cell_Size EM-Muscle Cells
for ii = 1:M_Cell_Size
    EM_muscles_spike_timing{ii,1} = [];
end

%Sparse each spike and sort into array according to EM-Muscle Cell ID
%order.
for ii = 1:spike_count
    %Only EM-Muscle spikes are processed
    if spk.cellid(ii) >= 256 && spk.cellid(ii) <= (M_Cell_Size+255)
        %sorted numerically, starting from the 1st EM-Muscle cell (id 256)
        EM_muscles_spike_timing{spk.cellid(ii)-255,1} = [EM_muscles_spike_timing{spk.cellid(ii)-255,1} spk.tspike(ii)];
    end
end

%Embed the EM_muscle spike trains into outputs, 
embedded_data = {};
stim_duration = 200;

%Loop Through All Stim (Single Channel) files
if target_num == 1
    window_duration = 300; %duration of neural activity affected by a stim, e.g., 200ms for target 1.
elseif target_num == 3
    window_duration = 200; %duration of neural activity affected by a stim, e.g., 200ms for target 3.
end

window_duration_forward = window_duration; %duration of neural activity affected by a stim, e.g., 400ms.
window_duration_backward = window_duration;

backward_window_start_time = max(0, stim_starting_time -  window_duration_backward);
for channel_i = 1:M_Cell_Size
    %each channel
    embedded_data.forward{1,channel_i} = -1*ones(Max_Spike_Count_per_Bin,window_duration_forward); %partitioned into 1ms intervals/bins, with 100 max spikes per interval/bin.
    embedded_data.backward{1,channel_i} = -1*ones(Max_Spike_Count_per_Bin,window_duration_forward);%partitioned into 1ms intervals/bins, with 100 max spikes per interval/bin.
    
    for jj = 1:window_duration_forward
        %forward corresponds to time interval 200ms+jj-1 : 600ms+jj-1 in channel
        %stimulation at jj index affects a window from 200 to 600, shifted by jj.
        forward_valid_spk_time_start = EM_muscles_spike_timing{channel_i,1}(EM_muscles_spike_timing{channel_i,1}>=stim_starting_time+jj-1); %find all causual spikes
        forward_valid_spk_time = forward_valid_spk_time_start(forward_valid_spk_time_start < stim_starting_time + window_duration_forward +jj-1);  %limit all causual spikes within a 400ms window
        
        
        spike_count = length(forward_valid_spk_time);
        if spike_count > 0
            embedded_data.forward{1,channel_i}(1:spike_count,jj) = (forward_valid_spk_time - (stim_starting_time+jj-1))'; %normalize the spike timing within each window
        end
        
        %Do the same for anticausual response of 400 ms window duration
        backward_valid_spk_time_start = EM_muscles_spike_timing{channel_i,1}(EM_muscles_spike_timing{channel_i,1}> backward_window_start_time+ jj-1);
        backward_valid_spk_time = backward_valid_spk_time_start(backward_valid_spk_time_start < stim_starting_time+jj-1);
        
        spike_count = length(backward_valid_spk_time);
        if spike_count > 0
            embedded_data.backward{1,channel_i}(1:spike_count,jj) = (backward_valid_spk_time - (jj-1))'; %normalize the spike timing within each window
        end
    end
end

%% Update codebook and error_book

y = zeros(192,stim_duration);

for n=1:stim_duration
    %training
    ii = 1:codebook_length;
    
    forward_distance = zeros(M_Cell_Size,codebook_length);
    backward_distance = zeros(M_Cell_Size,codebook_length);
    %Compute the spike distance between the input and the
    %codebook/dictionary
    
    for jj = 1:M_Cell_Size
        %compute the spike train distance matrix (M_Cell_Size x dictionary size)
        %spike_distance(Spike Matrix; Spike Vector, Time_Length,
        %Max_Spike_Count)
        forward_distance(jj,:) =  spike_distance(codebook.forward{1,jj},embedded_data.forward{1,jj}(:,n),400,100);
        backward_distance(jj,:) = spike_distance(codebook.backward{1,jj},embedded_data.backward{1,jj}(:,n),400,100);
    end
    
    %use sum kernel between M_Cell_Size channels
    average_forward_distance = sum(forward_distance)/M_Cell_Size;
    average_backward_distance = sum(backward_distance)/M_Cell_Size;
    
    y(:,n) = lr_qk*error_book(:,ii)*(exp(-k_size_f*average_forward_distance-k_size_b*average_backward_distance))';
end

[value, location] = max(y(:));
y_norm = y./value;
y_desired_norm = y_norm;
y_desired = y;

%% Generating firing rate patterns
% Compute discrete firing rates from the continuous signal, and convert to the format required by NEURON

incInterval = 10;  % min time interval to consider adding increasing values
decInterval = 50; % min time interval to consider adding decresing values
maxTime = 300; % max duration of stim
decThreshold = 0.3; % minimum decrease factor to consider adding value

spike_pattern = {};
spike_pattern_cord = [];
spike_pattern_Display_cord = [];
stim_i = 0;
for row_i = 1:length(y_norm(:,1))
    prev_rate = 0;
    prev_time = -300;
    cell = row_i -1;
    for col_i = 1:length(y_norm(1,:))
        rate = round(y_norm(row_i,col_i)*maxRate);
        time = col_i - 1 + 100;
        % if increased and min inc interval has passed, include value
        if rate > prev_rate && (time - prev_time) > incInterval
            stim_i = stim_i + 1;
            spike_pattern{1,stim_i}{1,1} = cell;
            spike_pattern{1,stim_i}{1,2} = time;
            spike_pattern{1,stim_i}{1,3} = maxTime;
            spike_pattern{1,stim_i}{1,4} = rate;
            spike_pattern_cord(stim_i,:) = [cell  time  rate];
            spike_pattern_Display_cord(stim_i,:) = [row_i,  col_i, y_norm(row_i,col_i)];
            prev_rate = rate;
            prev_time = time;
        % if decreased by minThreshold factor, and max inc interval has passed, include value
        elseif abs(rate - prev_rate) > (decThreshold*maxRate) && (time - prev_time) > decInterval
            stim_i = stim_i + 1;
            spike_pattern{1,stim_i}{1,1} = cell;
            spike_pattern{1,stim_i}{1,2} = time;
            spike_pattern{1,stim_i}{1,3} = maxTime;
            spike_pattern{1,stim_i}{1,4} = rate;
            spike_pattern_cord(stim_i,:) = [cell  time  rate];
            spike_pattern_Display_cord(stim_i,:) = [row_i,  col_i, y_norm(row_i,col_i)];
            prev_rate = rate;
            prev_time = time;
        end
    end
    % add last time
    if prev_time > 0
        stim_i = stim_i + 1;
        spike_pattern{1,stim_i}{1,1} = cell;
        spike_pattern{1,stim_i}{1,2} = maxTime;
        spike_pattern{1,stim_i}{1,3} = maxTime;
        spike_pattern{1,stim_i}{1,4} = 0;
        spike_pattern_cord(stim_i,:) = [cell  maxTime 0];
        spike_pattern_Display_cord(stim_i,:) = [row_i,  col_i, 0];
    end
end
 
%% Plot figure
plotFig=1;
if plotFig
    h1=figure
    [X_plot,Y_plot] = meshgrid(1:1:200);
    colors=-del2(250-(250*subplus(y_norm)));
    h = surf(X_plot(1:192,:),Y_plot(1:192,:),subplus(y_norm*250), colors);%, 'FaceLighting','gouraud');%,'LineWidth',0.3);%, 'facecolor', 'texturemap');

    title_info = ['Optimized Stim Estimate Across All Channels (Target ',num2str(target_num), ', Type ', num2str(perturbation_type), ', ',num2str(perturbation_percentage), '%)'];
    
    title(title_info);
    zlabel('Spike rate (Hz)', 'fontsize', 14);% x 250 Hz');
    ylabel('Cell', 'fontsize', 14);
    xlabel('Time (ms)', 'fontsize', 14);% (starting from 100 ms)');
    [value, location] = max(y_norm(:));
    [n_stim,n_col] = size(spike_pattern_cord);
    dx = 0;
    dy = 0;
    for stim_i = 1:n_stim
        b = ['Cell ID: ',num2str(spike_pattern_cord(stim_i,1)),', Start Time (ms): ',num2str(spike_pattern_cord(stim_i,2)), ', Rate (Hz): ', num2str(spike_pattern_cord(stim_i,3))]; c = cellstr(b);
    end
    %colorbar
    axis square
    set(gca,'XTickLabel',[100:50:300])
    zlim([0 250]);
    view([-37.5, 40])
    set(h1,'PaperUnits','inches','PaperPosition',[0 0 8 8])
    print(h1,'-dpng','-r600',[dataFolder dataSubfolder 'T',num2str(target_num),'T',num2str(perturbation_type),'P',num2str(perturbation_percentage) '_',probing,'.png']);
    
end
%% Write stim to txt file
fileID = fopen([dataFolder dataSubfolder 'target-',num2str(target_num),'_ptype-',num2str(perturbation_type),'_pperc-',num2str(perturbation_percentage)  '_' probing  '_repair-in.txt'],'w');
[n_stim,n_col] = size(spike_pattern_cord);
for stim_i = 1:n_stim
    %start time
    fprintf(fileID,'%2.3f %2.3f %2.3f\r\n',spike_pattern{1,stim_i}{1,1},spike_pattern{1,stim_i}{1,2},spike_pattern{1,stim_i}{1,4});
    %end time
end
fclose(fileID);
