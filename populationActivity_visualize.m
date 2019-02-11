%% read in sorted spike matfiles, get statistics from each recording 
%  and finally combine into master statistics per condition (WT vs. KO)


%% set up paths 

results_folder = uigetdir();

datFiles = dir(fullfile(results_folder,'*.mat'));
datFiles = {datFiles(:).name};

metadatFiles = dir(fullfile(results_folder,'*.txt'));
metadatFiles = {metadatFiles(:).name};

%% loop through data files and plot population histograms

all_histograms = {};

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
    
    % quick check on whether T is correctly computed -- if not, then use max timestamp recorded
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
        
    % create population histogram
    
    bin_width = 0.1; % bin width, in seconds
    
    hist_edges = 0:bin_width:T;
    
    pop_hist = smooth(histcounts(all_tmsp./Fs,hist_edges,'Normalization','count'),5);
    
    x_axis = hist_edges(1:end-1) + diff(hist_edges)./2;
        
    hist_array = [pop_hist,x_axis'];
    
    if ~isempty(strfind(mat_fnam,'Control')) || ~isempty(strfind(mat_fnam,'WT'))
        all_histograms = [all_histograms;{[1],mat_fnam, hist_array}];
    elseif ~isempty(strfind(mat_fnam,'Mutant')) || ~isempty(strfind(mat_fnam,'KO'))
        all_histograms = [all_histograms;{[2],mat_fnam, hist_array}];
    end
    
end

%% plot population histograms of all cultures' data


% first, sort histograms by condition (flag in the first column of
% all_histograms, 1 indicates WT, 2 indicates KO)

[~,srt_idx] = sort(cell2mat(all_histograms(:,1)),'ascend');
all_histograms = all_histograms(srt_idx,:);


% indices of histograms to plot
demonstrative_hists_idx = [1 3 4 6 7 9 10 11 12 13 15 16 17 19 20 21 22 24 25];
filtered_histograms = all_histograms(demonstrative_hists_idx,:);


% create spacing (y-axis) for histogram plots
num_hists_2plot = length(filtered_histograms);
inter_hist_height = 20; % distance between population histograms on figure;
height_vector = 1:inter_hist_height:(num_hists_2plot*inter_hist_height);


figure(2);
for f_i = 1:size(filtered_histograms,1)
    
    switch filtered_histograms{f_i,1} 
        case 1
            plot(filtered_histograms{f_i,3}(:,2),filtered_histograms{f_i,3}(:,1)+height_vector(f_i),'b-','LineWidth',1);
        case 2
            plot(filtered_histograms{f_i,3}(:,2),filtered_histograms{f_i,3}(:,1)+height_vector(f_i),'r-','LineWidth',1);
    end
    
    hold on;
    
end

axis tight;

xlabel('Time (seconds)','FontSize',15)
ylabel(sprintf('Activity (spikes / %d milliseconds)',bin_width * 1000),'FontSize',15)

title('Population activity for wild-type vs. mutant cortical cultures','FontSize',18)




    
    