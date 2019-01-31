%% set paths, using readNexFile to invoke file-choice UI

% 06-12-2018, Conor Heins

fprintf('Please choose working directory\n');
working_dir = uigetdir();
cd(working_dir);
addpath(fullfile(working_dir,'NexReader'))
addpath(fullfile(working_dir,'plotSpikeRaster_v1.2'));

fprintf('Please choose .nex file desired for analysis\n');
[fnam,fdir] = uigetfile('*.nex');

nexFile = readNexFile(fullfile(fdir,fnam));

%% read in data, using 'Filter' flag to selectively read in filtered data

% parameters for spike extraction
sigma_rms = 4; % multiples of standard deviation of the noise to use for spike detection
Fs = nexFile.freq; % sampling frequency (for Sindhu, 20 kHz)
pre_ms = 2; % pre-spike time in milliseconds
post_ms = 4; % post-spike time in milliseconds
dead_ms = 2; % dead-time between successive spike detections

filter_flag = true;
high_pass_freq = 100; % in Hertz

if filter_flag
    cont_var_idx = determine_channels(nexFile,'contvars','Electrode');
    
    % extract spike waveforms and exclude large voltage artifacts (with
    % optional inspection afterwards)
    tic
    [waveforms,timestamps] = extract_waveforms(nexFile,'contvars',cont_var_idx,sigma_rms,Fs,dead_ms,pre_ms,post_ms,filter_flag,high_pass_freq);
    fprintf('Time taken to extract waveforms from continuous trace: %.2f seconds\n',toc)

else
    cont_var_idx = determine_channels(nexFile,'contvars','Filter');
    
    % extract spike waveforms and exclude large voltage artifacts (with
    % optional inspection afterwards)
    tic
    [waveforms,timestamps] = extract_waveforms(nexFile,'contvars',cont_var_idx,sigma_rms,Fs,dead_ms,pre_ms,post_ms,filter_flag,high_pass_freq);
    fprintf('Time taken to extract waveforms from continuous trace: %.2f seconds\n',toc)
    
end

ignore_channels = find(cellfun(@isempty,waveforms));
    
waveforms(ignore_channels) = [];
timestamps(ignore_channels) = [];

%% artifact removal (based on exceeding standard deviations of average waveform, user-chosen SD-factor)

artifact_removal_flag = false;

if artifact_removal_flag
    
    SD_exclude_factor = 15;
    [waveforms,timestamps,num_wf,excluded_WF] = remove_artifacts(waveforms,timestamps,SD_exclude_factor);
    
    figure(1);
    display_waveforms(waveforms,cellfun(@(x) ~isempty(x), excluded_WF))
    close gcf;
end

%% cubic spline interpolation and peak-alignment

spline_alignment = false;

if spline_alignment
    tic
    dt = 1/Fs; % sampling interval of the raw data
    upsample_factor = 4; % factor by which to interpolate the waveform
    sample_idx = 36:44; % indices in time (original sampling) to interpolate -- ideally, one would choose indices that surround the peak
    
    [waveforms_aligned,waveforms_splined] = interpolate_and_align_sindhu(waveforms,dt,upsample_factor,sample_idx);
    fprintf('Time taken to interpolate & peak-align: %.2f seconds\n',toc);
    
    
    figure(2);
    display_waveforms(waveforms_aligned);
    close gcf;
    
    waveforms = waveforms_aligned;

end

tic
window_indices = { [32:38], [38:46], [50:80] }; % windows chosen to extract pre-spike polarization, peri-spike depolarization, and post-spike 'AHP'-thing
thresholds = [1e-3,3e-2]; % standard deviation and spike amplitude waveforms used to finally classify waveform as 
                         % a real spike or not
[spike_waveforms,spike_tmsp,noise_wf] = threshold_waveforms_sindhu(waveforms,timestamps,window_indices,thresholds);
fprintf('Time taken to classify spikes as true or false positives: %.2f seconds\n',toc)

%% plot spike rasters

% prepare data structure representation for plotting
active_channels = find(cellfun(@(x) ~isempty(x), spike_tmsp));
spike_tmsp_active = spike_tmsp(active_channels)';


for ii= 1:length(active_channels)
    spike_tmsp_active{ii} = (spike_tmsp_active{ii}/Fs)'; % convert to seconds by dividing by sampling frequency
end

[xPoints, yPoints] = plotSpikeRaster(spike_tmsp_active,'PlotType','vertline');

recording_nam = fnam(1:(end-4));
title(sprintf('Population activity from %s',recording_nam))
xlabel('Time (seconds)')
ylabel('Channel ID')
axis tight;

%% analysis

active_channels = find(cellfun(@(x) ~isempty(x), spike_tmsp));
spike_tmsp_active = spike_tmsp(active_channels)';

spk = [];
for chan = 1:length(active_channels)
    spk_tmp = spike_tmsp_active{chan}';
    spk = [spk; [repmat(active_channels(chan),size(spk_tmp,1),1), spk_tmp]];
end

[~,idx] = sort(spk(:,2),'ascend');
spk = spk(idx,:);

[population_activity,time_x] = hist(spk(:,2),1000);

% get firing rates, coefficients of variation, and correlation matrices

T = size(nexFile.contvars{1}.data,1)/Fs;
bin_size = 0.01; % bin-size, in seconds
bin_edges = 0:bin_size:T;

firing_rates = zeros(length(active_channels),1);
cvs = zeros(length(active_channels),1);
correlations = zeros(length(active_channels));

for unit1 = 1:length(active_channels)
    
    spk_temp = spk(spk(:,1) == active_channels(unit1),2);
    
    binned1 = histc(spk_temp,bin_edges);
    firing_rates(unit1) = length(spk_temp)/T;
    cvs(unit1) = std(diff(spk_temp))./mean(diff(spk_temp));
    
    
    for unit2 = 1:length(active_channels)
        
        binned2 = histc(spk(spk(:,1) == active_channels(unit2),2),bin_edges);
        
        correlations(unit1,unit2) = corr(binned1,binned2);
    end
    
end






    



 
