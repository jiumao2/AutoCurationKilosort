function [cluster_id_A, cluster_id_B] = splitCluster(folder_data, cluster_id, spike_ids_A)
% SPLITCLUSTER Split a cluster into two clusters as phy does.
%
% Input:
%   - folder_data: the folder where the data is located
%   - cluster_id: 1x1 value. The cluster id to process
%   - spike_ids_A: nx1 vector. The spikes to be removed. These spikes should
%   be in the given cluster.
%
% Output:
%   - cluster_id_A: 1x1 value. 
%   - cluster_id_B: 1x1 value.
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
spike_ids_all = find(spike_clusters == cluster_id);

% check all the spike_ids are in the cluster
assert(length(intersect(spike_ids_all, spike_ids_A)) == length(spike_ids_A));
assert(~isempty(spike_ids_A));
assert(length(spike_ids_all) > length(spike_ids_A));

% modify the spike_clusters
spike_ids_B = setdiff(spike_ids_all, spike_ids_A);

max_cluster_id = max(spike_clusters);
cluster_id_A = max_cluster_id+1;
cluster_id_B = max_cluster_id+2;
spike_clusters(spike_ids_A) = cluster_id_A;
spike_clusters(spike_ids_B) = cluster_id_B;

% deal with cluster_groups
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    idx_this = find(cluster_group.cluster_id == cluster_id);
    if ~isempty(idx_this)
        group_this = cluster_group.group{idx_this};
        cluster_group(idx_this,:) = [];
        
        tbl = table();
        tbl.cluster_id = [cluster_group.cluster_id; cluster_id_A; cluster_id_B];
        tbl.group = [cluster_group.group; group_this; group_this];
        writetable(tbl, fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
    end
end

% deal with KSLabels
cluster_KSLabel = readtable(fullfile(folder_data, 'cluster_KSLabel.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_this = find(cluster_KSLabel.cluster_id == cluster_id);
if ~isempty(idx_this)
    KSLabel_this = cluster_KSLabel.KSLabel{idx_this};
    cluster_KSLabel(idx_this,:) = [];
    tbl = table();
    tbl.cluster_id = [cluster_KSLabel.cluster_id; cluster_id_A; cluster_id_B];
    tbl.KSLabel = [cluster_KSLabel.KSLabel; KSLabel_this; KSLabel_this];
    writetable(tbl, fullfile(folder_data, 'cluster_KSLabel.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end

% deal with cluster_Amplitude
cluster_Amplitude = readtable(fullfile(folder_data, 'cluster_Amplitude.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_this = find(cluster_Amplitude.cluster_id == cluster_id);
if ~isempty(idx_this)
    Amplitude_this = cluster_Amplitude.Amplitude(idx_this);
    cluster_Amplitude(idx_this,:) = [];
    
    tbl = table();
    tbl.cluster_id = [cluster_Amplitude.cluster_id; cluster_id_A; cluster_id_B];
    tbl.Amplitude = [cluster_Amplitude.Amplitude; Amplitude_this; Amplitude_this];
    writetable(tbl, fullfile(folder_data, 'cluster_Amplitude.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end

% deal with cluster_ContamPct
cluster_ContamPct = readtable(fullfile(folder_data, 'cluster_ContamPct.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_this = find(cluster_ContamPct.cluster_id == cluster_id);
if ~isempty(idx_this)
    ContamPct_this = cluster_ContamPct.ContamPct(idx_this);
    cluster_ContamPct(idx_this,:) = [];
    
    tbl = table();
    tbl.cluster_id = [cluster_ContamPct.cluster_id; cluster_id_A; cluster_id_B];
    tbl.ContamPct = [cluster_ContamPct.ContamPct; ContamPct_this; ContamPct_this];
    writetable(tbl, fullfile(folder_data, 'cluster_ContamPct.tsv'), 'Delimiter', '\t', 'FileType', 'text');
end

% save the spike_clusters
writeNPY(spike_clusters, fullfile(folder_data, 'spike_clusters.npy'));

end