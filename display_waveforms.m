function [ ] = display_waveforms( waveforms,subset_idx )
%display_waveforms 
% Quick wrapper for some plotting functions for displaying
% candidate spike waveforms stored in {num chans x 1} cell array
%   Can provide subset_idx to only display certain waveforms

if ~exist('subset_idx','var') || isempty(subset_idx)
    
    for i = 1:length(waveforms)
        plot(waveforms{i}');
        pause;
    end
    
else
    
    for i = 1:length(waveforms)
        plot(waveforms{i}(subset_idx,:)');
        pause;
    end 
    
end

end

