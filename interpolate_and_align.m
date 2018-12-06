function [ waveforms_aligned,waveforms_splined,proc_time ] = interpolate_and_align(waveforms,SA_params)
%interpolate_and_align_sindhu 
% Upsampling and cubic spline interpolation to align peaks of spike
% waveforms within given temporal range of the segment (provided by 'sample_idx')
%   Need to provide {numchans x 1} cell array of waveforms, sampling
%   interval dt, upsample_factor and indices of where to perform spline
%   interpolation within

tic

try, sample_idx = SA_params.sample_idx; catch, sample_idx = 36:44; end
try, dt= SA_params.dt; catch, dt = 1/20e3; end
try, upsample_factor = SA_params.upsample_factor; catch, upsample_factor = 4; end
try, plot_flag = SA_params.plot_flag; catch, plot_flag = true; end


numchans = length(waveforms);

xt = [0 : (length(sample_idx) - 1) ] .* dt; % time domain of original data
xxt = linspace(0,max(xt), length(xt) * upsample_factor); % time domain of interpolation;

% interpolation

waveforms_splined = cell(numchans,1);

for chan = 1:numchans
    temp_mat = waveforms{chan};
    
    if size(temp_mat,2) == 1
        temp_mat = temp_mat';
    end
   
    num_wf = size(temp_mat,1);
    temp_mat_spline = zeros(num_wf,length(sample_idx)*upsample_factor);
    for wf = 1:num_wf
        temp_wf = temp_mat(wf,sample_idx);
        temp_wf_spline = spline(xt,temp_wf,xxt);
        temp_mat_spline(wf,:) = temp_wf_spline;
    end
    waveforms_splined{chan} = temp_mat_spline;
    
end

% alignment

waveforms_aligned = cell(numchans,1);

for chan = 1:numchans
    temp_mat_spline = waveforms_splined{chan};
    temp_mat_aligned = zeros(size(waveforms{chan}));
    temp_mat_spline(:,[1:30, 50:end]) = 0; % make sure to look only around peak -- this re-designed for Sindhu's data (20 kHz)
    [~,peak_idx] = min(temp_mat_spline,[],2); %indices of the peaks
    idx_to_align = mode(peak_idx);
    offsets = idx_to_align - peak_idx;

    for wf = 1:numel(peak_idx)
        [~,original_idx] = min(abs(xt - xxt(peak_idx(wf))));
        [~,shifted_idx] = min(abs(xt - xxt(peak_idx(wf) + offsets(wf))));
        shift_wf = shifted_idx - original_idx;
        if shift_wf ~= 0
            temp_mat_aligned(wf,:) = circshift(waveforms{chan}(wf,:),shift_wf,2);
        else
            temp_mat_aligned(wf,:) = waveforms{chan}(wf,:);
        end
    end
    
    waveforms_aligned{chan} = temp_mat_aligned;
end

if plot_flag
    figure;
    display_waveforms(waveforms_aligned);
    pause;
end

proc_time = toc;


end

