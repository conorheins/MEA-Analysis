
% good example of wild-type data -- f_i is the index of the datFiles cell
% array
% f_i = 19

% hist(all_tmsp,1000);

bin_width = 0.2; % bin width, in seconds

hist_edges = 0:bin_width:T;

% pop_hist = smooth(histcounts(all_tmsp./Fs,hist_edges,'Normalization','count'),5);
pop_hist = histcounts(all_tmsp./Fs,hist_edges,'Normalization','count');

x_axis = hist_edges(1:end-1) + diff(hist_edges)./2;

bar(x_axis,pop_hist,1.5,'k')
xlim([0 600])

ax = gca;
ax.FontSize = 16;
xlabel('Time (seconds)','FontSize',16)
ylabel('Population activity (spikes per 200 ms bin)','FontSize',16)



burst_colors = jet(size(NB,1));
for i = 1:size(NB,1)
    
    hold on; 
    
    plot([NB(i,1) NB(i,1)],[0 100],'--','LineWidth',1.5,'Color',burst_colors(i,:));
    plot([NB(i,2) NB(i,2)],[0 100],'--','LineWidth',1.5,'Color',burst_colors(i,:));
    
end