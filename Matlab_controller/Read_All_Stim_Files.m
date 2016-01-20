close all
clear all
clc

load('iterparams.mat');

%% params
M_Cell_Size = 192;
stim_duration = 200;
rate_multiplier = 1;
starting_time_multiplier = 1;

% read the existing inverse model codebook and errorbook
load([dataFolder dataSubfolder 'Inv_Plant_target-',num2str(target_num),'_ptype-',num2str(perturbation_type),'_pperc-',num2str(perturbation_percentage)  '_' probing '.mat']);

%Loop Through All Stim (Single Channel) files
if target_num == 1
    window_duration = 300; %duration of neural activity affected by a stim, e.g., 200ms for target 1.
elseif target_num == 3
    window_duration = 200; %duration of neural activity affected by a stim, e.g., 200ms for target 3.
end

window_duration_forward = window_duration; %duration of neural activity affected by a stim, e.g., 400ms.
window_duration_backward = window_duration;

for cell_id = 1:M_Cell_Size-1
    cell_id
    codebook_length
    rate = rate_multiplier*250; %spiking frequency
    stim_start_time = starting_time_multiplier*100; %stim starting time 200ms or 400ms
    
    %Read from this data file:
    file_name = ['target-',num2str(target_num),'_ptype-',num2str(perturbation_type),'_pperc-',num2str(perturbation_percentage),'_cell-',num2str(cell_id-1),'_start-',num2str(stim_start_time),'_dur-',num2str(stim_duration),'_rate-',num2str(rate)];
    file_type = '_spk.mat';
    
    dname = bmmFolderPerturb;
    fname=fullfile(dname,[file_name,file_type]);
    spk = load(fname);
    
    %% Record the spiking activity in the M_Cell_Size Muscle-driving EM population:
    EM_muscles_spike_timing = {};
    %spike_count = length(spk.tspike); %total number of spikes in the set, i.e., total number of spikes across 663 cells in 1.02s duration
    target_reaching_spikes = spk.tspike;
    if target_num == 1
        target_reaching_spikes = target_reaching_spikes(target_reaching_spikes <= 600);
    elseif target_num == 3
        target_reaching_spikes = target_reaching_spikes(target_reaching_spikes <= 500);
    end
    spike_count = length(target_reaching_spikes);
    
    %Allocation spike timing array for each of the M_Cell_Size EM-Muscle Cells
    for ii = 1:M_Cell_Size
        EM_muscles_spike_timing{ii,1} = [];
    end
    
    %Sparse each spike and sort into array according to EM-Muscle Cell ID order.
    for ii = 1:spike_count
        %Only EM-Muscle spikes are processed
        if spk.cellid(ii) >= 256 && spk.cellid(ii) <= (M_Cell_Size+255)
            %sorted numerically, starting from the 1st EM-Muscle cell (id 256)
            EM_muscles_spike_timing{spk.cellid(ii)-255,1} = [EM_muscles_spike_timing{spk.cellid(ii)-255,1} spk.tspike(ii)];
        end
    end
    
    %%
    %Data Formatting for Output (EM-Muscle) Response
    %========================================================================================================================
    %Embed the EM_muscle spike trains into outputs
    embedded_data = {};
    
    backward_window_start_time = max(0, stim_start_time -  window_duration_backward);
    for channel_i = 1:M_Cell_Size
        %each channel
        embedded_data.forward{1,channel_i} = -1*ones(Max_Spike_Count_per_Bin,window_duration_forward); %partitioned into 1ms intervals/bins, with 100 max spikes per interval/bin.
        embedded_data.backward{1,channel_i} = -1*ones(Max_Spike_Count_per_Bin,window_duration_forward);%partitioned into 1ms intervals/bins, with 100 max spikes per interval/bin.
        
        %    length = EM_muscles_spike_timing{ii,1}
        %    for EM_muscles_spike_timing{ii,1}  = 1
        %
        %    end
        for jj = 1:window_duration_forward
            %forward corresponds to time interval 200ms+jj-1 : 600ms+jj-1 in channel
            %stimulation at jj index affects a window from 200 to 600, shifted
            %by jj.
            forward_valid_spk_time_start = EM_muscles_spike_timing{channel_i,1}(EM_muscles_spike_timing{channel_i,1}>=stim_start_time+jj-1); %find all causual spikes
            forward_valid_spk_time = forward_valid_spk_time_start(forward_valid_spk_time_start < stim_start_time + window_duration_forward +jj-1);  %limit all causual spikes within a 400ms window
            
            
            spike_count = length(forward_valid_spk_time);
            if spike_count > 0
                embedded_data.forward{1,channel_i}(1:spike_count,jj) = (forward_valid_spk_time - (stim_start_time+jj-1))'; %normalize the spike timing within each window
            end
            
            %Do the same for anticausual response of 400 ms window duration
            backward_valid_spk_time_start = EM_muscles_spike_timing{channel_i,1}(EM_muscles_spike_timing{channel_i,1}> backward_window_start_time+ jj-1);
            backward_valid_spk_time = backward_valid_spk_time_start(backward_valid_spk_time_start < stim_start_time+jj-1);
            
            
            spike_count = length(backward_valid_spk_time);
            if spike_count > 0
                embedded_data.backward{1,channel_i}(1:spike_count,jj) = (backward_valid_spk_time - (jj-1))'; %normalize the spike timing within each window
            end
        end
    end
    %%  Update codebook and error_book
    
    %Desired Signal
    desired_stim = zeros(M_Cell_Size,stim_duration);
    desired_stim(cell_id, (stim_start_time - 100 +1):stim_start_time + 100) = (rate/rate)*ones(1,200);
    
    y = zeros(M_Cell_Size,stim_duration);
    %start
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
            forward_distance(jj,:) = spike_distance(codebook.forward{1,jj},embedded_data.forward{1,jj}(:,n),400,100);
            backward_distance(jj,:) = spike_distance(codebook.backward{1,jj},embedded_data.backward{1,jj}(:,n),400,100);
        end
        
        %use sum kernel between M_Cell_Size channels
        average_forward_distance = sum(forward_distance)/M_Cell_Size;
        average_backward_distance = sum(backward_distance)/M_Cell_Size;
        
        y(:,n) = lr_qk*error_book(:,ii)*(exp(-k_size_f*average_forward_distance-k_size_b*average_backward_distance))';
        error = desired_stim(:,n) - y(:,n);
        
        %Quantization
        %M_Cell_Size channels and the max distance from all channels is recorded
        forward_min = max(sqrt(forward_distance));
        backward_min = max(sqrt(backward_distance));
        
        quantize = 0;
        Index_min = 0;
        for jj = 1:codebook_length
            if forward_min(jj) <= q_factor && backward_min(jj) <= q_factor
                quantize = 1;
                Index_min = jj;
                break
            end
        end
        
        if quantize == 1
            error_book(:,Index_min) = error_book(:,Index_min) + error;
        else
            for jj = 1:M_Cell_Size
                codebook.forward{1,jj} = [codebook.forward{1,jj} embedded_data.forward{1,jj}(:,n)];
                codebook.backward{1,jj} = [codebook.backward{1,jj} embedded_data.backward{1,jj}(:,n)];
            end
            codebook_length = codebook_length+1;
            error_book = [error_book zeros(M_Cell_Size,1)];
            error_book(:,n) = error;
        end
    end
end

save([dataFolder dataSubfolder  'Inv_Plant_target-',num2str(target_num),'_ptype-',num2str(perturbation_type),'_pperc-',num2str(perturbation_percentage) '_' probing '.mat'],'error_book', 'codebook', 'codebook_length', 'q_factor', 'k_size_f', 'k_size_b', 'lr_qk','Max_Spike_Count_per_Bin','target_num','perturbation_type', 'perturbation_percentage');
