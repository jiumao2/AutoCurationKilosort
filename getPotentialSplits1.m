function getPotentialSplits(folder_data, user_settings)
% GETPOTENTIALSPLITS find potential splits
%
% Input:
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
%
% Output:
%   To do
%
% Determination of potential splits:
%   (1) Exclude all noise clusters
%   (2) Any sub-cluster should contain at least 10% of spikes
%   (3) Not continuous in time and PETH is different in the two segments
%   (4) Two peaks are found in the amplitude-time figure, which can be
%   separated using a timg-varing gaussian model.
%   (5) Amplitude in the waveform are separatable and could be separated
%   using a time-varing gaussian model.
%

% read the params
min_cluster_percentage = user_settings.splitting.min_cluster_percentage;

verbose = user_settings.merging.verbose;

% load the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));

cluster_ids = unique(spike_clusters);

% get the non-noise clusters
cluster_id_noise = [];
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    cluster_id_noise = cluster_group.cluster_id(strcmpi(cluster_group.group, 'noise'));
end

cluster_id_non_noise = setdiff(cluster_ids, cluster_id_noise);


%% dectect potential splits with dbscan
split_info = {}; % cluster_id, spike_ids_A, spike_ids_B

for k = 1:length(cluster_id_non_noise)
    spike_ids = find(spike_clusters == cluster_id_non_noise(k));
    amplitudes_this = amplitudes(spike_ids);
    spike_times_this_sec = double(spike_times(spike_ids))./30000;

    min_num_spikes = 5*length(spike_times_this_sec)/100;

    mdl = TimeVaryingGaussian();
    mdl.fit(amplitudes_this, spike_times_this_sec);
    p = mdl.pdf(amplitudes_this, spike_times_this_sec);
    p_crit = normpdf(3);

    idx = p < p_crit;
    if sum(idx==1) < max(min_num_spikes, 500)
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
        continue
    end

    spike_ids_A = spike_ids(idx==0);
    spike_ids_B = spike_ids(idx==1);
    
    split_info{end+1} = {cluster_id_non_noise(k), spike_ids_A, spike_ids_B};

    if verbose
        title_name = ['Cluster#' num2str(cluster_id_non_noise(k)),...
            ' (n = ', num2str(length(spike_times_this_sec)), ')'];
        save_filename = fullfile(folder_data, 'Fig/PotentialSplits/',...
            ['Cluster' num2str(cluster_id_non_noise(k))]);
        
        plotClusterComparison(folder_data,...
            spike_ids_A,...
            spike_ids_B,...
            title_name,...
            save_filename);
    end
    
    fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
end



