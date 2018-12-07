%% Batch Processing for Pre-Processing and Spike-Sorting MEA data stored in .nex files
% 06-12-2018, Conor Heins
% last edit, 07-12-2018

%% set up paths and directories
fprintf('Please choose working directory\n');
working_dir = uigetdir();
cd(working_dir);
addpath(fullfile(working_dir,'NexReader'))
addpath(fullfile(working_dir,'plotSpikeRaster_v1.2'));

%% choose folder to write to

fprintf('Please choose results master folder to write to\n')
results_master_dir = uigetdir();

%% set up parameters for processing

% parameters for pre-processing
processParams.field_nam = 'contvars'; % choose field within nexFile struct to get e-phys data from
processParams.sigma_rms = 4; % multiples of standard deviation of the noise to use for spike detection
processParams.pre_ms = 2; % pre-spike time in milliseconds
processParams.post_ms = 4; % post-spike time in milliseconds
processParams.dead_ms = 2; % dead-time between successive spike detections
processParams.filter_flag = true;   % when importing data, using 'Filter' flag to selectively process the already-filtered data
if processParams.filter_flag
    processParams.high_pass_freq = 100; % in Hertz
end

artifact_removal_flag = false; % set to true if you want to perform automatic artifact detection/removal
if artifact_removal_flag
    AR_params.SD_exclude_factor = 15; % factor of the standard deviation to use for artifact detect detection
    AR_params.plot_flag = true;
end

spline_alignment = false; % set to true if you want to use cubic interpolation for better alignment of spike waveforms
if spline_alignment
    SA_params.upsample_factor = 4;  % factor by which to interpolate the waveform
    SA_params.sample_idx = 36:44; % indices in time (original sampling) to interpolate -- ideally, one would choose indices that surround the peak
    SA_params.plot_flag = true; 
end

sortParams.window_indices = { [32:38], [38:46], [50:80] }; % windows chosen to extract pre-spike polarization, peri-spike depolarization, and post-spike 'AHP'-thing
sortParams.thresholds = [1e-3,3e-2]; % standard deviation and spike amplitude thresholds for waveforms (in mV), used to finally classify waveform as a real spike or not

%% load data for analysis

[master_dir,filenames] = choose_files();

text_gen_function = @(string_i,number_i)(sprintf('%s\t%.2f seconds',string_i,number_i));
proc_stages = {'Data import';'Waveform/Tmsp extraction';'Artifact Removal';'Interp./Peak-Alignment';'Spike Sorting'};

results_fdir = fullfile(results_master_dir,sprintf('results_%s_%s',datestr(now,'mmddyy'),datestr(now,'HH')));

if exist(results_fdir,'dir') ~= 7
    mkdir(results_fdir)
end

%% loop through files
for f_i = 1:length(filenames)
    
    processing_times = cell(5,1);
    
    %% first step, read in the data (.nex file read function)
    tic
    nexFile= readNexFile(fullfile(master_dir,filenames{f_i}));
    time_taken = toc;
    fprintf('Time taken to read .nex file into workspace: %.2f seconds\n',time_taken)
    processing_times{1} = time_taken;
    processParams.Fs = nexFile.freq;
    
    %% second step, extract waveforms from electrode traces
    [wfs,tmsp,names,time_taken] = extract_waveforms(nexFile,processParams);
    fprintf('Time taken to extract waveforms from continuous trace: %.2f seconds\n',time_taken)
    processing_times{2} = time_taken;
    
    duration_Seconds = length(nexFile.contvars{1}.data)/processParams.Fs;
    clear nexFile;
    
    % get rid of channels in which no waveforms were detected
    ignore_channels = find(cellfun(@isempty,wfs));
    
    wfs(ignore_channels) = [];
    tmsp(ignore_channels) = [];
    names(ignore_channels) = [];
    
    
    %% third step, perform optional artifact removal 

    if artifact_removal_flag
        [wfs,tmsp,num_wf,excluded_WF,time_taken] = remove_artifacts(wfs,tmsp,AR_params);
        processing_times{3} = time_taken;  
        fprintf('Time taken to remove artifacts: %.2f seconds\n',time_taken);
    else
        processing_times{3} = 0;
    end
        
    %% fourth step, perform optional cubic spline interpolation and peak-alignment

    if spline_alignment        
        SA_params.dt = 1/processParams.Fs; % sampling interval of the raw data
        [wfs,waveforms_splined,time_taken] = interpolate_and_align(wfs,SA_params);
        processing_times{4} = time_taken;
        fprintf('Time taken to interpolate & peak-align: %.2f seconds\n',time_taken);
    else
        processing_times{4} = 0;
    end
    
    %% fifth step, sort spikes
    
    [spike_waveforms,spike_tmsp,noise_wf,time_taken] = threshold_waveforms(wfs,tmsp,sortParams);
    fprintf('Time taken to classify spikes as true or false positives: %.2f seconds\n',time_taken)
    processing_times{5} = time_taken;
    
    
    %% save sorted spikes results, create metadata file, etc.

    non_empty_idx = find(~cellfun(@isempty,spike_waveforms));
    for chan_i = 1:length(non_empty_idx)
        spikes(chan_i).filename = fullfile(master_dir,filenames{f_i});        
       
        if exist('extractBetween.m','file') ~= 2
            split_up = strsplit(names{non_empty_idx(chan_i)},'_');
            number_idx = find(~isnan(cellfun(@str2double,split_up)));
            spikes(chan_i).Label = str2double(split_up{number_idx(1)});
            if number_idx > 1
                spikes(chan_i).ID = str2double(split_up{number_idx(2)});
            else
                spikes(chan_i).ID = 'Ref';
            end
        else
            spikes(chan_i).Label = str2double(extractBetween(names{non_empty_idx(chan_i)},'Label_','_ID'));
            spikes(chan_i).ID = str2double(extractBetween(names{non_empty_idx(chan_i)},'ID_','_Str'));
        end
        
        spikes(chan_i).waveforms = spike_waveforms{non_empty_idx(chan_i)};
        spikes(chan_i).timestamps = spike_tmsp{non_empty_idx(chan_i)};
    end
    
    components = strsplit(filenames{f_i},filesep);
    textFID = fopen(fullfile(results_fdir,sprintf('metadata_%s.txt',components{end}(1:end-4))),'w');
    fprintf(textFID,'Full path to nex file: %s\n',fullfile(master_dir,filenames{f_i}));
    lines2write = cellfun(@(str_i,num_i)text_gen_function(str_i,num_i),proc_stages,processing_times,'UniformOutput',false);
    fprintf(textFID,'%s\n',lines2write{:});
    fprintf(textFID,'Number of active electrodes: %d\n',length(non_empty_idx));
    fprintf(textFID,'Length of recording: %.2f seconds\n',duration_Seconds);
    fclose(textFID);
    
    all_vars = who;
    param_var_indices = find(~cellfun(@isempty,strfind(all_vars,'Params')));
    param_vars = all_vars(param_var_indices);
    
    vars_to_save = ['spikes';param_vars];
    
    save(fullfile(results_fdir,sprintf('%s_sorted.mat',components{end}(1:end-4))),vars_to_save{:},'-v7.3');
    
    clearvars spikes;
         
      
end
  
   
    



    
    
