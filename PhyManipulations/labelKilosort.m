function labelKilosort(folder_data, cluster_ids, labels)
% LABELKILOSORT Labels the clusters as 'good', 'mua' or 'noise' as phy does.
%
% Input:
%   - folder_data: the folder where the data is located
%   - cluster_ids: nx1 vector. The cluster IDs which will be updated.
%   - labels: char or nx1 cell. The labels of the clusters: 'good', 'mua', 'noise' or 'unsorted'.
%       The 'unsorted' labels will not saved in cluster_group.tsv.
%
% Output:
%   The cluster_group.tsv will be updated.
%

% check the input
if size(cluster_ids, 1) == 1
    cluster_ids = cluster_ids';
end

if ~iscell(labels)
    label_this = labels;
    labels = cell(length(cluster_ids), 1);
    for k = 1:length(cluster_ids)
        labels{k} = label_this;
    end
end

if size(labels, 1) == 1
    labels = labels';
end

% read the old cluster groups
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
cluster_id_old = [];
group_old = {};
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    cluster_id_old = cluster_group.cluster_id;
    group_old = cluster_group.group;
end

% remove the cluster_id that will be updated
idx_kept = arrayfun(@(id) ~any(cluster_ids == id), cluster_id_old);
cluster_id_old = cluster_id_old(idx_kept);
group_old = group_old(idx_kept);

% combined to create new cluster ids and groups
cluster_id_new = [cluster_id_old; cluster_ids];
group_new = [group_old; labels];

% remove the unsorted labels
idx_unsorted = find(strcmpi(group_new, 'unsorted'));
cluster_id_new(idx_unsorted) = [];
group_new(idx_unsorted) = [];

% save to cluster_group.tsv
tbl = table();
tbl.cluster_id = cluster_id_new;
tbl.group = group_new;

writetable(tbl, fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end