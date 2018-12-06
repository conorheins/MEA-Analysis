function [waveforms,timestamps,num_wf,artifacts,proc_time] = remove_artifacts(waveforms,timestamps,AR_params)
% remove_artifacts
%   removes artifacts based on traces that exceed SD_exclude_factor * standard deviation of
%   the mean voltage value at any temporal sample

tic;
try SD_exclude_factor = AR_params.SD_exclude_factor; catch, SD_exclude_factor = 10;end
try plot_flag = AR_params.plot_flag; catch, plot_flag = true; end

numchans = length(timestamps);
num_wf = zeros(numchans,1);
artifacts = cell(numchans,1);

for chan = 1:numchans
    temp_wf = waveforms{chan};
    std_dev_wf = std(temp_wf,0,1);
    mean_wf = mean(temp_wf,1);
    outlier_thr = [mean_wf - SD_exclude_factor * std_dev_wf; mean_wf + SD_exclude_factor*std_dev_wf];
    [~,bad_ids] = find(or(temp_wf < outlier_thr(1), temp_wf > outlier_thr(2)));
    bad_ids = unique(bad_ids);
    artifacts{chan} = waveforms{chan}(bad_ids,:);
    waveforms{chan}(bad_ids,:) = [];
    timestamps{chan}(bad_ids) = [];
    num_wf(chan) = size(waveforms{chan},1);
end

if plot_flag   
    figure;
    display_waveforms(waveforms,cellfun(@(x) ~isempty(x), artifacts))
    pause;
    close gcf;
end

finish = toc;

proc_time = finish;

end

