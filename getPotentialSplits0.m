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
n_categories = user_settings.splitting.n_categories;
dbscan_minpts = user_settings.splitting.dbscan_minpts;
dbscan_epsilon = ...
    user_settings.splitting.dbscan_epsilon_range(1):user_settings.splitting.dbscan_epsilon_step:user_settings.splitting.dbscan_epsilon_range(2);
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
[minpts_mesh, epsilon_mesh] = meshgrid(dbscan_minpts, dbscan_epsilon);

n_mdl = length(minpts_mesh(:));

for k = 1:length(cluster_id_non_noise)
    spike_ids = find(spike_clusters == cluster_id_non_noise(k));
    amplitudes_this = amplitudes(spike_ids);
    spike_times_this = double(spike_times(spike_ids));

    min_num_spikes = min_cluster_percentage*length(spike_times_this)/100;

    data = [zscore(spike_times_this), zscore(amplitudes_this)];
    
    silhouette_values = NaN(1, n_mdl);
    idx_all = cell(1, n_mdl);
    for j = 1:n_mdl
        idx = dbscan(data, epsilon_mesh(j), minpts_mesh(j));
        if any(max(idx) == n_categories)   
            idx_unlabeled = find(idx<0);
            idx_labeled = find(idx>0);

            if sum(idx==1) < min_num_spikes || sum(idx==2) < min_num_spikes
                continue
            end
            
            mdl = fitcknn(data(idx_labeled,:), idx(idx_labeled), 'NumNeighbors', minpts_mesh(j));
            temp = idx;
            temp(idx_unlabeled) = mdl.predict(data(idx_unlabeled,:));
            silhouette_values(j) = median(silhouette(data, temp));
            idx_all{j} = {spike_ids(temp==1), spike_ids(temp==2)};
        end
    end

    if all(isnan(silhouette_values))
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
        continue 
    end

    [~, idx_max] = max(silhouette_values);
    spike_ids_A = idx_all{idx_max}{1};
    spike_ids_B = idx_all{idx_max}{2};
    
    split_info{end+1} = {cluster_id_non_noise(k), spike_ids_A, spike_ids_B};

    if verbose
        title_name = ['Cluster#' num2str(cluster_id_non_noise(k)),...
            ' (n = ', num2str(length(spike_times_this)), ')'];
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



