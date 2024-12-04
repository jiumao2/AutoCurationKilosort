function updateClusterInfo(folder_data)
% UPDATECLUSTERINFO Auto update cluster_info.tsv as phy does
% 
% Input:
%   - folder_data: the folder where the data is located
% 
% Output:
%   updated cluster_info.tsv
% 
% Variables in cluster_info.tsv:
% (1) cluster_id: from phy
% (2) Amplitude: from Kilosort
% (3) ContamPct: from Kilosort
% (4) KSLabel: from Kilosort
% (5) amp: the median computed from amplitudes.npy
% (6) ch: the channel with the largest amplitude computed from the templates
% (7) depth: the y coordinate of the channel
% (8) fr: mean firing rate
% (9) group: manually labeled
% (10) n_spikes: number of spikes
% (11) sh: filled with 0
% 

disp('Updating cluster_info.tsv...');

% read all npy files
amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));
% channel_map = readNPY(fullfile(folder_data, 'channel_map.npy'));
channel_positions = readNPY(fullfile(folder_data, 'channel_positions.npy'));
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
spike_templates = readNPY(fullfile(folder_data, 'spike_templates.npy'));
templates = readNPY(fullfile(folder_data, 'templates.npy'));
templates_ind = readNPY(fullfile(folder_data, 'templates_ind.npy'));

% read all tsv files
cluster_Amplitude = table2array(readtable(fullfile(folder_data, 'cluster_Amplitude.tsv'), 'Delimiter', '\t', 'FileType', 'text'));
cluster_ContamPct = table2array(readtable(fullfile(folder_data, 'cluster_ContamPct.tsv'), 'Delimiter', '\t', 'FileType', 'text'));
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
cluster_KSLabel = readtable(fullfile(folder_data, 'cluster_KSLabel.tsv'), 'Delimiter', '\t', 'FileType', 'text');

%% compute the needed data
cluster_id = unique(spike_clusters);
duration_sec = double((max(spike_times) - min(spike_times)))./30000;
n_spike = arrayfun(@(id)sum(spike_clusters == id), cluster_id);
fr = n_spike./duration_sec;

n_templates = size(templates, 1);
[~, idx_max] = max(squeeze(max(templates, [], 2) - min(templates, [], 2)), [], 2);
ch_templates_ind0 = arrayfun(@(k) templates_ind(k, idx_max(k)), (1:n_templates)');
cluster_templates_ind0 = arrayfun(@(id)mode(spike_templates(spike_clusters == id)), cluster_id);
ch_ind0 = ch_templates_ind0(cluster_templates_ind0+1);
depth = channel_positions(ch_ind0+1, 2);

% group
group = cell(length(cluster_id), 1);
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    for k = 1:length(cluster_id)
        idx = find(cluster_group.cluster_id == cluster_id(k));
        if isempty(idx)
            continue
        end

        group{k} = cluster_group.group{idx};
    end
end

tbl = table();
tbl.cluster_id = cluster_id;

% Amplitude
idx_amp = findSeq(cluster_Amplitude(:,1), cluster_id);
tbl.Amplitude = cluster_Amplitude(idx_amp, 2);

% ContamPct
idx_contam = findSeq(cluster_ContamPct(:,1), cluster_id);
tbl.ContamPct = cluster_ContamPct(idx_contam, 2);

% KSLabel
idx_kslabel = findSeq(cluster_KSLabel.cluster_id, cluster_id);
tbl.KSLabel = cluster_KSLabel.KSLabel(idx_kslabel);

tbl.amp = arrayfun(@(id)median(amplitudes(spike_clusters == id)), cluster_id);
tbl.ch = ch_ind0;
tbl.depth = depth;
tbl.fr = fr;
tbl.group = group;
tbl.n_spike = n_spike;
tbl.sh = zeros(length(cluster_id), 1);

writetable(tbl, fullfile(folder_data, 'cluster_info.tsv'), 'Delimiter', '\t', 'FileType', 'text');

disp('Cluster_info.tsv is updated!');
end