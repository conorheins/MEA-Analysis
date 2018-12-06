function [ bad_ids ] = artifact_remove_singlechannel(wfs_single_channel,SD_factor)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

outlier_thr = SD_factor * std(wfs_single_channel,0,1); % assumes individual waveforms are in the rows, temporal dimensions/samples in the columns
[~,bad_ids] = find(or(wfs_single_channel < mean(wfs_single_channel,1) - max(outlier_thr), wfs_single_channel > mean(wfs_single_channel,1) + max(outlier_thr)));

end

