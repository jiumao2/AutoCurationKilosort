t_scale = 1;

st = double(spike_times(spike_clusters == 663 | spike_clusters == 1095));
amp = double(amplitudes(spike_clusters == 663 | spike_clusters == 1095));

st_zscore = normalize(st, 'zscore');
amp_zscore = normalize(amp, 'zscore');

%%
data = double([st_zscore*t_scale, amp_zscore]);
distanceMatrix = pdist(data);

k_all = [];
idx_all = {};
silhouette_value_all = [];
for k = 1:-0.05:0.01
    idx = dbscan(data, k, 50);
    if max(idx) == 2
        k_all = [k_all, k];

        idx_unlabeled = find(idx<0);
        idx_labeled = find(idx>0);
        
        mdl = fitcknn(data(idx_labeled,:), idx(idx_labeled), 'NumNeighbors', 50);
        temp = idx;
        temp(idx_unlabeled) = mdl.predict(data(idx_unlabeled,:));
        silhouette_value_all = [silhouette_value_all, median(silhouette(data, temp))];

        idx_all{end+1} = idx;
    end
end

fig = EasyPlot.figure();
ax_all = EasyPlot.createGridAxes(fig, 1, length(k_all));

for k = 1:length(ax_all)
    plot(ax_all{k}, data(idx_all{k}==1, 1), data(idx_all{k}==1, 2), 'b.');
    plot(ax_all{k}, data(idx_all{k}==2, 1), data(idx_all{k}==2, 2), 'r.');
    plot(ax_all{k}, data(idx_all{k}<=0, 1), data(idx_all{k}<=0, 2), 'k.');
    title(ax_all{k}, ['k = ',num2str(k_all(k)), ' (', num2str(silhouette_value_all(k)), ')']);
end

EasyPlot.cropFigure(fig);