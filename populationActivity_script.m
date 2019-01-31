%% read in sorted spike matfiles, get statistics from each recording 
%  and finally combine into master statistics per condition (WT vs. KO)


%% set up paths 

addpath(genpath('burst_netBurst_detection_algo'));

results_folder = uigetdir();

datFiles = dir(fullfile(results_folder,'*.mat'));
datFiles = {datFiles(:).name};

metadatFiles = dir(fullfile(results_folder,'*.txt'));
metadatFiles = {metadatFiles(:).name};

%% loop through data files and analyze spike-trains and bursts


WT_SingleChannels = {}; % data array to eventually be exported into csvs/excel spreadsheets
KO_SingleChannels = {}; % arrays to eventually be exported into csvs/excel spreadsheets

WT_NetworkBursts = {}; % data array to eventually be exported into csvs/excel spreadsheets
KO_NetworkBursts = {}; % arrays to eventually be exported into csvs/excel spreadsheets

for f_i = 1:length(datFiles)
    
    mat_fnam = datFiles{f_i};
    load(fullfile(results_folder,mat_fnam));
    
    fnamPatt = strsplit(mat_fnam,'_');
    fnamPatt = fnamPatt{1};
    
    metadat_fnam = fullfile(results_folder,metadatFiles{contains(metadatFiles,fnamPatt)});
    
    fileID = fopen(fullfile(results_folder,metadatFiles{contains(metadatFiles,fnamPatt)}),'r');
    
    blah = '';
    while ~strcmp(blah,'recording:')
        blah = fscanf(fileID,'%s\t.2f');
    end
    
    duration_str = fscanf(fileID,'%s\t.2f');
    fclose(fileID);
    T = str2double(duration_str);
    
    Fs = processParams.Fs;
    
    numChannels = length(spikes);
    
    %quick check on whether T is correctly computed -- if not, then use max timestamp recorded
    % (This was done because of weird metadata file for recording:
    % 2018-06-08T15-52-2408062018, that indicated recording only lasted 0.8
    % seconds, whereas spike timestamps indicate at least a 15 minute
    % recording
    
    all_tmsp = [];
    for chan_i = 1:numChannels
        all_tmsp = [all_tmsp;spikes(chan_i).timestamps];
    end
    
    if T*Fs < max(all_tmsp)
        T = max(all_tmsp)./Fs;
    end
        
    
    %% single channel analysis 
    
    channel_stats = zeros(numChannels,4); 
    % per channel statistics (1st column: firing rates; 2nd column: coefficient
    % of variation of ISI distribution; 3rd column: burst rate; 4th column: coefficient of variation
    % of IBI distribution
    
    all_spikes = [];
    all_bursts_sparse_vec = spalloc(T*Fs,1,T*Fs);
    BDTrain = [];
    
    burst_stats = zeros(3,4); % burst statistics (1st column: mean/std/n of network burst rate; 2nd column: mean/std/n of 
    % total number bursts per network burst; 3rd column: mean/std/n of
    % network burst duration; 4th column: mean/std/n of total number
    % participating channels per network burst)
    
    for chan_i = 1:numChannels
        
        spike_times = unique(spikes(chan_i).timestamps);
        
        if length(spike_times) >= 5
            all_spikes = [all_spikes;[repmat(spikes(chan_i).ID,length(spike_times),1), spike_times./Fs]];
            
            %1. calculate firing rate as # number of spikes / total duration of
            %recording
            channel_stats(chan_i,1) = length(spike_times)./T;
            
            %2. calculate coefficient of variation of ISI distribution
            
            ISIs = diff([0;spike_times./Fs;T]);
            channel_stats(chan_i,2) = std(ISIs)./mean(ISIs);
            
            %% burst detection algorithm (based on Pasquale et al. 2010, logISI method)
            
            % put spike_times into a big sparse vector
            
            spike_times_sparse = sparse(spike_times,ones(length(spike_times),1),ones(length(spike_times),1),T*Fs,1,T*Fs);
            
            nBinsPerDec = 10;
            [bins, ISIlog_hist_norm, allISI] = calcISILogHist(spike_times_sparse, nBinsPerDec, Fs); % Valentina's code
            
            mpd = 2;
            voidParamTh = 0.7;
            ISITh = 150;
            [ISImax, ~, flags] = calcISImax(bins, ISIlog_hist_norm, mpd, voidParamTh, ISITh); % Valentina's code
            
            if ~flags(2)
                ISImax = ISITh; % if calcISImax fails to find well-separated peaks, then default to ISITh
            end
            
            ISImaxTh = ISImax + ceil(ISImax/5);
            maxISImaxTh = ISImax + ceil(ISImax/5);
            % NOTES ABOUT THESE TWO PARAMETERS ^
            % if Condition 1. ISImax < maxISImaxTh AND Condition 2. ISImax > ISImaxTh are TRUE, then ISImaxTh is used as the intra-burst spike detector,
            % and ISImax is used to extend the intra-burst window to include boundary spikes. Else, if only Condition 1 is satisfied, then ISImax is
            % used as the intra-burst spike detector and there is no boundary spike inclusion. And if NOT even Condition 1 is satisfied, then
            % ISImaxTh is used as the intra-burst detector, regardless of whether it's smaller or larger than ISImax. And in this case, also no boundary spike
            % inclusion.
            minNumSpikes = 3;
            [burstTrain, burstEventTrain] = burstDetection(spike_times_sparse,ISImax,flags, ISImaxTh, maxISImaxTh, Fs, minNumSpikes);
            
            if ~isempty(burstTrain)
                all_IBIs = burstTrain(:,5);
                channel_stats(chan_i,3) = 1./mean(all_IBIs);
                channel_stats(chan_i,4) = std(all_IBIs)./mean(all_IBIs);
                
                % accumulate single-channel burst information into array-wide burst
                % information
                
                all_bursts_sparse_vec(find(burstEventTrain),1) = 1;
                
                nbursts_temp = size(burstTrain,1);
                BDTrain = [BDTrain; [repmat(spikes(chan_i).ID,nbursts_temp,1), burstTrain(:,1:2)] ];
            else
                channel_stats(chan_i,3:4) = NaN;
            end
        else
            channel_stats(chan_i,:) = NaN;
        end
      
    end
    
    % now do network burst detection, using accumulated single-channel
    % bursts
    
    if numChannels < 5
        
        if ~isempty(strfind(mat_fnam,'Control')) || ~isempty(strfind(mat_fnam,'WT'))
            WT_SingleChannels = [WT_SingleChannels;{mat_fnam, channel_stats}];          
        elseif ~isempty(strfind(mat_fnam,'Mutant')) || ~isempty(strfind(mat_fnam,'KO'))
            KO_SingleChannels = [KO_SingleChannels;{mat_fnam, channel_stats}];            
        end
        
    else
    
        nBinsPerDec = 10;
        [bins, IBeIhist, allIBI] = calcISILogHist(all_bursts_sparse_vec, nBinsPerDec, Fs); % Valentina's code
        
        mpd = 2;
        voidParamTh = 0.7;
        IBeITh = 1000;
        [IBeImax, pks, flags] = calcIBeImax(bins, IBeIhist, mpd, voidParamTh, IBeITh); % Valentina's code
        
        IBeIThDef = 250;
        
        if ~flags(2)
            IBeImax = IBeIThDef;
        end
        
        [NB] = netBurstDetection(BDTrain, IBeImax, flags, IBeIThDef, Fs, floor(numChannels*0.25)); % Valentina's code
        
        %% Fill out network burst statistics
        
        % 1. network burst rates
        nb_rates = 1./diff(NB(:,1)./Fs);
        burst_stats(1,1) = mean(nb_rates);
        burst_stats(2,1) = std(nb_rates);
        burst_stats(3,1) = length(nb_rates);
        
        % 2. single channel bursts / network burst
        burst_stats(1,2) = mean(NB(:,3));
        burst_stats(2,2) = std(NB(:,3));
        burst_stats(3,2) = length(NB(:,3));
        
        % 3. network burst duration
        burst_stats(1,3) = mean(NB(:,4)./Fs);
        burst_stats(2,3) = std(NB(:,4)./Fs);
        burst_stats(3,3) = length(NB(:,4));
        
        % 4. number of participating electrodes / burst
        burst_stats(1,4) = mean(NB(:,5));
        burst_stats(2,4) = std(NB(:,5));
        burst_stats(3,4) = length(NB(:,5));
        
        
        %% assemble statistics into master arrays that accumulate data across conditions
        
        if ~isempty(strfind(mat_fnam,'Control')) || ~isempty(strfind(mat_fnam,'WT'))
            WT_SingleChannels = [WT_SingleChannels;{mat_fnam, channel_stats}];
            WT_NetworkBursts = [WT_NetworkBursts;{mat_fnam,burst_stats}];
            
        elseif ~isempty(strfind(mat_fnam,'Mutant')) || ~isempty(strfind(mat_fnam,'KO'))
            KO_SingleChannels = [KO_SingleChannels;{mat_fnam, channel_stats}];
            KO_NetworkBursts = [KO_NetworkBursts;{mat_fnam,burst_stats}];
            
        end
    end
        
end

%% accumulate all data into summary arrays to copy and paste in Excel -- SINGLE CHANNEL SUMMARIES


%% Wild-Types
recording_names = cell(size(WT_SingleChannels,1),1);

single_channel_summary = zeros(size(WT_SingleChannels,1),8);

for i = 1:size(WT_SingleChannels,1)
    
    recording_names{i} = WT_SingleChannels{i,1}(1:end-11);
    
    single_channel_summary(i,1) = size(WT_SingleChannels{i,2},1); % Num Active Channels
    
    FRs = WT_SingleChannels{i,2}(:,1);
    FRs = FRs(~isnan(FRs));
    single_channel_summary(i,2) = mean(FRs); % Average Firing Rate
    single_channel_summary(i,3) = median(FRs); % Median Firing Rate
    
    CVs = WT_SingleChannels{i,2}(:,2);
    CVs = CVs(~isnan(CVs));
    single_channel_summary(i,4) = mean(CVs); % Mean Coefficient of Variation
    single_channel_summary(i,5) = median(CVs); % Median Coefficient of Variation
    
    IBIs = WT_SingleChannels{i,2}(:,2);
    IBIs = IBIs(~isnan(IBIs));
    single_channel_summary(i,6) = mean(IBIs); % Mean Inter-Burst Interval
    single_channel_summary(i,7) = median(IBIs); % Median Inter-Burst Interval
    single_channel_summary(i,8) = std(IBIs)/mean(IBIs); % Coefficient of IBI Distribution
    
end

%% Knock-Outs
recording_names = cell(size(KO_SingleChannels,1),1);

single_channel_summary2 = zeros(size(KO_SingleChannels,1),8);

for i = 1:size(KO_SingleChannels,1)
    
    recording_names{i} = KO_SingleChannels{i,1}(1:end-11);
    
    single_channel_summary2(i,1) = size(KO_SingleChannels{i,2},1); % Num Active Channels
    
    FRs = KO_SingleChannels{i,2}(:,1);
    FRs = FRs(~isnan(FRs));
    single_channel_summary2(i,2) = mean(FRs); % Average Firing Rate
    single_channel_summary2(i,3) = median(FRs); % Median Firing Rate
    
    CVs = KO_SingleChannels{i,2}(:,2);
    CVs = CVs(~isnan(CVs));
    single_channel_summary2(i,4) = mean(CVs); % Mean Coefficient of Variation
    single_channel_summary2(i,5) = median(CVs); % Median Coefficient of Variation
    
    IBIs = KO_SingleChannels{i,2}(:,2);
    IBIs = IBIs(~isnan(IBIs));
    single_channel_summary2(i,6) = mean(IBIs); % Mean Inter-Burst Interval
    single_channel_summary2(i,7) = median(IBIs); % Median Inter-Burst Interval
    single_channel_summary2(i,8) = std(IBIs)/mean(IBIs); % Coefficient of IBI Distribution
    
end

%% accumulate all data into summary arrays to copy and paste in Excel -- NETWORK BURST SUMMARIES


%% Wild-Types
recording_names = cell(size(WT_NetworkBursts,1),1);

NB_summary = zeros(size(WT_NetworkBursts,1),3);

for i = 1:size(WT_NetworkBursts,1)
    
    recording_names{i} = WT_NetworkBursts{i,1}(1:end-11);
    
    num_bursts = WT_NetworkBursts{i,2}(3,:);
    num_bursts = num_bursts(~isnan(num_bursts));
    NB_summary(i,1) = max(num_bursts);
    
    NB_summary(i,2) = WT_NetworkBursts{i,2}(1,1);
    NB_summary(i,3) = WT_NetworkBursts{i,2}(1,2);
    NB_summary(i,4) = WT_NetworkBursts{i,2}(1,3);
    NB_summary(i,5) = WT_NetworkBursts{i,2}(1,4);
   
    
end

%% Knock-Outs

recording_names = cell(size(KO_NetworkBursts,1),1);

NB_summary2 = zeros(size(KO_NetworkBursts,1),3);

for i = 1:size(KO_NetworkBursts,1)
    
    recording_names{i} = KO_NetworkBursts{i,1}(1:end-11);
    
    num_bursts = KO_NetworkBursts{i,2}(3,:);
    num_bursts = num_bursts(~isnan(num_bursts));
    NB_summary2(i,1) = max(num_bursts);
    
    NB_summary2(i,2) = KO_NetworkBursts{i,2}(1,1);
    NB_summary2(i,3) = KO_NetworkBursts{i,2}(1,2);
    NB_summary2(i,4) = KO_NetworkBursts{i,2}(1,3);
    NB_summary2(i,5) = KO_NetworkBursts{i,2}(1,4);
   
    
end
    

%%
    
    
    




   
           
           
           
           
           
           
           
           
           
           
            
            
            
            
            
            
            
            
            
        
        
        
        
        
        
        
        
        
        
  
    
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

