function cluster_id_new = mergeCluster(folder_data, cluster_ids)
% MERGECLUSTER Merge clusters into one cluster as phy does.
%
% Input:
%   - folder_data: the folder where the data is located
%   - cluster_ids: nx1 value. The cluster ids to merge
%
% Output:
%   - cluster_id_new: 1x1 value. 
%
% The files will be changed:
% (1) spike_clusters.npy
% (2) cluster_group.tsv
% (3) cluster_Amplitude.tsv
% (4) cluster_ContamPct.tsv
% (5) cluster_KSLabel.tsv
%

% read spike_clusters
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));

% modify the spike_clusters
max_cluster_id = max(spike_clusters);
cluster_id_new = max_cluster_id+1;
for k = 1:length(cluster_ids)
    spike_clusters(spike_clusters == cluster_ids(k)) = cluster_id_new;
end

% deal with cluster_groups
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    idx_this = [];
    for k = 1:length(cluster_ids)
        temp = find(cluster_group.cluster_id == cluster_ids(k));
        if ~isempty(temp)
            idx_this = [idx_this; temp];
        end
    end

    if ~isempty(idx_this)
        group_this = cluster_group.group(idx_this);
        if any(strcmpi(group_this, 'noise'))
            group_out = 'noise';
        elseif any(strcmpi(group_this, 'mua'))
            group_out = 'mua';
        else
            group_out = 'good';
        end

        cluster_group(idx_this,:) = [];
        
        tbl = table();
        tbl.cluster_id = [cluster_group.cluster_id; cluster_id_new];
        tbl.group = [cluster_group.group; group_out];
        writetable(tbl, fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
    end
end

% deal with KSLabels
cluster_KSLabel = readtable(fullfile(folder_data, 'cluster_KSLabel.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_this = [];
for k = 1:length(cluster_ids)
    temp = find(cluster_KSLabel.cluster_id == cluster_ids(k));
    if ~isempty(temp)
        idx_this = [idx_this; temp];
    end
end

if ~isempty(idx_this)
    KSLabel_this = cluster_KSLabel.KSLabel(idx_this);
    if any(strcmpi(KSLabel_this, 'noise'))
        KSLabel_out = 'noise';
    elseif any(strcmpi(KSLabel_this, 'mua'))
        KSLabel_out = 'mua';
    else
        KSLabel_out = 'good';
    end

    cluster_KSLabel(idx_this,:) = [];
    tbl = table();
    tbl.cluster_id = [cluster_KSLabel.cluster_id; cluster_id_new];
    tbl.KSLabel = [cluster_KSLabel.KSLabel; KSLabel_out];
    writetable(tbl, fullfile(folder_data, 'cluster_KSLabel.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end

% deal with cluster_Amplitude
cluster_Amplitude = readtable(fullfile(folder_data, 'cluster_Amplitude.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_this = [];
for k = 1:length(cluster_ids)
    temp = find(cluster_Amplitude.cluster_id == cluster_ids(k));
    if ~isempty(temp)
        idx_this = [idx_this; temp];
    end
end

if ~isempty(idx_this)
    Amplitude_this = cluster_Amplitude.Amplitude(idx_this);
    cluster_Amplitude(idx_this,:) = [];
    
    tbl = table();
    tbl.cluster_id = [cluster_Amplitude.cluster_id; cluster_id_new];
    tbl.Amplitude = [cluster_Amplitude.Amplitude; mode(Amplitude_this)];
    writetable(tbl, fullfile(folder_data, 'cluster_Amplitude.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end

% deal with cluster_ContamPct
cluster_ContamPct = readtable(fullfile(folder_data, 'cluster_ContamPct.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_this = [];
for k = 1:length(cluster_ids)
    temp = find(cluster_ContamPct.cluster_id == cluster_ids(k));
    if ~isempty(temp)
        idx_this = [idx_this; temp];
    end
end

if ~isempty(idx_this)
    ContamPct_this = cluster_ContamPct.ContamPct(idx_this);
    cluster_ContamPct(idx_this,:) = [];
    
    tbl = table();
    tbl.cluster_id = [cluster_ContamPct.cluster_id; cluster_id_new];
    tbl.ContamPct = [cluster_ContamPct.ContamPct; mode(ContamPct_this)];
    writetable(tbl, fullfile(folder_data, 'cluster_ContamPct.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end

% save the spike_clusters
writeNPY(spike_clusters, fullfile(folder_data, 'spike_clusters.npy'));

end