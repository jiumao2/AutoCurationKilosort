function [cluster_id_new, cluster_id_removed] = removeSpikes(folder_data, cluster_id, spike_ids)
% REMOVESPIKES Remove spikes in a cluster.
%
% Input:
%   - folder_data: the folder where the data is located
%   - cluster_id: 1x1 value. The cluster id to process
%   - spike_ids: nx1 vector. The spikes to be removed. These spikes should
%   be in the given cluster.
%
% Output:
%   - cluster_id_removed: 1x1 value. The removed spikes generates a new
%   cluster.
%   - cluster_id_new: 1x1 value.
% 
% The files will be changed:
% (1) spike_clusters.npy
% (2) cluster_group.tsv
% (3) cluster_Amplitude.tsv
% (4) cluster_ContamPct.tsv
% (5) cluster_KSLabel.tsv
%

[cluster_id_removed, cluster_id_new] = splitCluster(folder_data, cluster_id, spike_ids);
labelKilosort(folder_data, cluster_id_removed, 'noise');

end