function [sorted_wfs,sorted_tmsp,noise_wfs,time_taken] = threshold_waveforms(waveforms,timestamps,sortParams)

tic

try window_indices = sortParams.window_indices; catch, window_indices = { [32:38], [38:46], [50:80] }; end
try thresholds = sortParams.thresholds; catch, thresholds = [1e-3,3e-2]; end

num_chans = length(waveforms);
sorted_wfs = cell(1,num_chans);
sorted_tmsp = cell(1,num_chans);
noise_wfs = cell(1,num_chans);

pre_spike_idx = window_indices{1};
during_spike_idx = window_indices{2};
post_spike_idx = window_indices{3};

for chan = 1:num_chans
    
    temp_chan_data = waveforms{chan};
    noise_IDX = true(size(temp_chan_data,1),1);
    if size(temp_chan_data,1) > 5
      
        pre_spike_vals = prctile(temp_chan_data(:,pre_spike_idx),90,2);
        [pre_ind_clust,C] = kmeans(pre_spike_vals,2);
        [~,higher_ind] = max(C);
        good_pre = pre_ind_clust == higher_ind;
        
        during_spike_vals = prctile(temp_chan_data(:,during_spike_idx),10,2);
        [during_spike_clust,C] = kmeans(during_spike_vals,2);
        [~,higher_ind] = min(C);
        good_during = during_spike_clust == higher_ind;
        
        post_spike_vals = prctile(temp_chan_data(:,post_spike_idx),80,2);
        [post_spike_clust,C] = kmeans(post_spike_vals,2);
        [~,higher_ind] = max(C);
        good_post = post_spike_clust == higher_ind;
        
        real_spk_idx = good_during & or(good_pre,good_post); % must have good spike amplitude and one of the other two features also has to be good
        
        avg_spk_wf = mean(temp_chan_data(real_spk_idx,:),1);
        
        if and(std(avg_spk_wf) > thresholds(1), abs(min(avg_spk_wf)) > thresholds(2)) % standard deviation of mean wf must be greater than threshold(1), max greater than threshold(2)
            sorted_wfs{chan} = temp_chan_data(real_spk_idx,:);
            sorted_tmsp{chan} = timestamps{chan}(real_spk_idx);
            noise_IDX(real_spk_idx) = false;
        end
    end
    
    noise_wfs{chan} = [find(noise_IDX),temp_chan_data(noise_IDX,:)];
      
end

time_taken = toc;