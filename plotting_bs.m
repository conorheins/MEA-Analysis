
% good example of wild-type data -- f_i is the index of the datFiles cell
% array

% f_i = 19;

% good example of KO data -- f_i is the index of the datFiles cell
% array

f_i = 8;

results_folder = uigetdir();

datFiles = dir(fullfile(results_folder,'*.mat'));
datFiles = {datFiles(:).name};

metadatFiles = dir(fullfile(results_folder,'*.txt'));
metadatFiles = {metadatFiles(:).name};

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

figure(2);
bin_width = 0.2; % bin width, in seconds

hist_edges = 0:bin_width:T;

pop_hist = histcounts(all_tmsp./Fs,hist_edges,'Normalization','count');

x_axis = hist_edges(1:end-1) + diff(hist_edges)./2;

bar(x_axis,pop_hist,2,'r')
xlim([0 T])

hold on;
spike_height = 1;
% start_height = max(pop_hist)+1;
start_height = 50;
raster_heights = start_height:spike_height:(start_height + (numChannels*spike_height));


% colors for WTs
colors = flipud(gray(2*numChannels));
colors = colors(7:(numChannels+6),:);

% colors for KOs
% colors = flipud(hot(3*numChannels));
% colors = colors(20:(numChannels+19),:);

for u_i = 1:numChannels
    
    [xOut, yOut] = rasterize(spikes(u_i).timestamps./Fs, raster_heights(u_i), raster_heights(u_i+1));
    
    plot(xOut,yOut,'Color',colors(u_i,:));
    
    hold on;
end
ylim([0 max(raster_heights)]);


ax = gca;
ax.FontSize = 16;
xlabel('Time (seconds)','FontSize',16)
ylabel('Population activity (spikes per 200 ms bin)','FontSize',16)
yticklabels = get(ax,'YTickLabel');
yticklabels{end} = 'Raster Plot'; % needs to exist but make it empty
set(ax,'YTickLabel',yticklabels)

title('Population Activity for Exemplary KO Culture','FontSize',18);