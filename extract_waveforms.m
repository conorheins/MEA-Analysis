function [ wfs,tmsp,names,proc_time,online_detect_flag ] = extract_waveforms(nexFile,processParams)
%extract_waveforms: From voltage time-series, use threshold-crossing and
%user-given  pre and post millisecond windows to pull out extracellular waveforms
% and put into two cell arrays: wfs (containing the waveforms from each
% channel) and tmsp (containing the timestamps from each channel)

tic

try, field_name = processParams.field_nam; catch, field_name = 'contvars'; end
try, Fs = processParams.Fs; catch, Fs = nexFile.freq; end
try, dead_ms = processParams.dead_ms; catch, dead_ms = 2; end
try, pre_ms = processParams.pre_ms; catch, pre_ms = 2; end
try, post_ms = processParams.post_ms; catch, post_ms = 2; end
try, sigma_rms = processParams.sigma_rms; catch, sigma_rms = 4; end
try, filter_flag = processParams.filter_flag; catch, filter_flag = true; end


if ~isfield(nexFile,field_name)
   field_name = 'waves';
   [channel_idx,names] = determine_channels(nexFile,field_name,'Spike Detector');
   wfs = cell(length(channel_idx),1);
   tmsp = cell(length(channel_idx),1);
   
   all_chans = nexFile.(field_name);
   for chan_id = 1:length(channel_idx)    
        wfs{chan_id} = all_chans{chan_id}.waveforms';
        tmsp{chan_id} = all_chans{chan_id}.timestamps * Fs;        
   end
   
   online_detect_flag = true;
   
   finish = toc;
   proc_time = finish;
   
   return
   
end

if filter_flag
    try, high_pass_fs = processParams.high_pass_freq; catch, high_pass_fs = 100; end
    [channel_idx,names] = determine_channels(nexFile,field_name,'Electrode');
else
    [channel_idx,names] = determine_channels(nexFile,field_name,'Filter');
end

dead_time = (dead_ms/1000) * Fs; %convert from millseconds (unit of dead_time) to samples
pre_time = (pre_ms/1000) * Fs;
post_time = (post_ms/1000) * Fs;

all_chans = nexFile.(field_name);

wfs = cell(length(channel_idx),1);
tmsp = cell(length(channel_idx),1);

for chan_id = 1:length(channel_idx)
    
    data = all_chans{channel_idx(chan_id)}.data;
    
    if filter_flag
        Wn = high_pass_fs / (Fs/2); % normalize high-pass cut-off frequency
        [b,a] = butter(2,Wn,'high');
        data = filtfilt(b,a,data);
    end
    
    thr = mean(data) - sigma_rms*std(data,1);
    
    crossings = find(data < thr); % less than threshold, because voltage crossings are negative in extracellular potentials
    
    if ~isempty(crossings)
        inter_cross_i = diff([0;crossings;length(data)]);
        
        excludeIDs = find(inter_cross_i <= dead_time);
        
        if excludeIDs(end) == (length(crossings)+1)
            crossings(excludeIDs(1:end-1)) = [];
            crossings(end) = [];
        else
            crossings(excludeIDs) = [];
        end
        
        num_candidates = length(crossings);
        
        temp_wfs = zeros(num_candidates,(pre_time + post_time + 1));
        
        throw_outs = false(num_candidates,1);
        
        for wf = 1:num_candidates
            if or(crossings(wf) - pre_time < 1, crossings(wf) + post_time > length(data))
                throw_outs(wf) = true;
            else
                segment = data( (crossings(wf) - pre_time) : (crossings(wf) + post_time) ); % cut out data segment
                [~,min_ind] = min(segment); % find minimum voltage within this segment
                crossings(wf) = crossings(wf) - pre_time + min_ind; % re-assignment timestamp to period of maximum
                if or(crossings(wf) - pre_time <= 0, crossings(wf) + post_time > length(data))
                    throw_outs(wf) = true;
                else
                    temp_wfs(wf,:) = data( (crossings(wf) - pre_time) : (crossings(wf) + post_time) );
                end
            end
        end
        
        temp_wfs(throw_outs,:) = [];
        crossings(throw_outs,:) = [];
        
        wfs{chan_id} = temp_wfs;
        tmsp{chan_id} = crossings;
    end
           
end

online_detect_flag =  false;


finish = toc;
proc_time = finish;



