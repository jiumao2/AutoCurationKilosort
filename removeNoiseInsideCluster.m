function removeNoiseInsideCluster(folder_data, user_settings)
% REMOVENOISEINSIDECLUSTER remove the "noise" in the cluster
%
% Determination of noise: outliers in the feature space. Utilize MATLAB
% isoutlier() function and PC features from Kilosort output to determine noise.
%
% Input (required):
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
% 
% Input (name-value pairs):
%   - verbose: false (default)
%       true / false. Control whether to generate figure output.
%   - n_pc_feature: 2 (default). 
%       Number of PC features to use. Range from 1 to 3.
%   - n_random_spikes: 100 (default)
%       The number of spikes to load from temp_wh.dat to compute the mean of the waveforms 
%       and determine the channel with the largest amplitude. 
%       It will be slow if the number is high.
%   - threshold: 5 (default)
%       The threshold used to determine the outlier in function isoutlier().
%       The code is `isoutlier(pc_features, 'median', 'ThresholdFactor', threshold)`;
%
% Output:
%   The outlier spikes will be spilted as new clusters and labeled as 'noise'. 
%   The corresponding files will be updated as phy does.
%

% parameters
n_pc_feature = user_settings.removeNoiseInsideCluster.n_pc_feature;
n_random_spikes = user_settings.removeNoiseInsideCluster.n_random_spikes;
threshold = user_settings.removeNoiseInsideCluster.threshold;
verbose = user_settings.removeNoiseInsideCluster.verbose;
waveform_window = user_settings.removeNoiseInsideCluster.waveform_window;

path_data = fullfile(folder_data, 'temp_wh.dat');
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_templates = readNPY(fullfile(folder_data, 'spike_templates.npy'));
pc_features = readNPY(fullfile(folder_data, 'pc_features.npy'));
pc_feature_ind = readNPY(fullfile(folder_data, 'pc_feature_ind.npy'));
load(fullfile(folder_data, 'ops.mat'));

% Get the IDs of non-noise clusters
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    cluster_noise = cluster_group.cluster_id(strcmpi(cluster_group.group, 'noise'));
else
    cluster_noise = [];
end

cluster_ids = unique(spike_clusters);
cluster_non_noise = setdiff(cluster_ids, cluster_noise);

% Clean the waveforms in non-noise clusters
disp('Start detecting the noise inside clusters!');
dir_output = dir(path_data);
nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});

for k = 1:length(cluster_non_noise)
    id = cluster_non_noise(k);
    spike_ids = find(spike_clusters == id);
    spike_times_this = spike_times(spike_ids);

    n_waveforms = min(length(spike_times_this), n_random_spikes);
    idx_rand = randperm(length(spike_times_this), n_waveforms);
    waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
    for j = 1:n_waveforms
        waveforms(j,:,:) = mmap.Data.x(:,...
            spike_times_this(idx_rand(j)) + waveform_window(1):spike_times_this(idx_rand(j)) + waveform_window(2));
    end
    
    mean_waveforms = squeeze(mean(waveforms, 1)); % 383 x 64
    [~, ch_largest] = max(max(mean_waveforms,[],2) - min(mean_waveforms,[],2));
    ch_largest_ind0 = ch_largest-1;

    % detect outliers
    is_outliers = false(1, length(spike_times_this));

    % load the pc features of the ch_largest
    pc_features_cluster = pc_features(spike_ids, :, :);
    pc_features_this = zeros(length(spike_times_this), n_pc_feature);
    spike_templates_this = spike_templates(spike_ids);
    for j = 1:length(spike_times_this)
        template_id_ind0 = spike_templates_this(j);
        ch = pc_feature_ind(template_id_ind0+1, :);
        idx = find(ch == ch_largest_ind0);
        if isempty(idx)
            is_outliers(j) = 1;
            pc_features_this(j,:) = NaN;
            continue
        end
        
        pc_features_this(j,:) = squeeze(pc_features_cluster(j, 1:n_pc_feature, idx));
    end
    
    for j = 1:n_pc_feature
        idx_found = isoutlier(pc_features_this(:,j), 'median', 'ThresholdFactor', threshold);
        is_outliers(idx_found) = 1;
    end

    if all(is_outliers == 0) || all(is_outliers == 1)
        continue
    end

    % plot the outliers
    if verbose
        % extract the waveforms of outliers
        n_outliers = sum(is_outliers);
        idx_outlier = find(is_outliers);
        waveforms_outliers = zeros(n_outliers, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
        for j = 1:n_outliers
            waveforms_outliers(j,:,:) = mmap.Data.x(:, spike_times_this(idx_outlier(j)) + waveform_window(1):spike_times_this(idx_outlier(j)) + waveform_window(2));
        end

        fig = EasyPlot.figure();
        ax_waveform = EasyPlot.axes(fig,...
            'Width', 6,...
            'Height', 6,...
            'XAxisVisible', 'off',...
            'YAxisVisible', 'off');
        ax_feature = EasyPlot.createAxesAgainstAxes(fig, ax_waveform, 'right',...
            'MarginBottom', 1,...
            'MarginLeft', 1);


        plot(ax_waveform, 1:diff(waveform_window)+1, squeeze(waveforms(:,ch_largest,:)), 'k-');
        plot(ax_waveform, 1:diff(waveform_window)+1, squeeze(waveforms_outliers(:,ch_largest,:)), 'r-');
        title(ax_waveform, ['Outliers: ', num2str(n_outliers), '/', num2str(length(spike_times_this))]);
        
        plot(ax_feature, pc_features_this(:,1), pc_features_this(:,2), 'k.');
        plot(ax_feature, pc_features_this(idx_outlier,1), pc_features_this(idx_outlier,2), 'r.');
        xlabel(ax_feature, 'PC1');
        ylabel(ax_feature, 'PC2');
        EasyPlot.cropFigure(fig);
        title(ax_feature, ['Cluster ', num2str(id)]);

        output_folder = fullfile(folder_data, 'Fig/RemoveWaveforms/');
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        output_filename = ['Cluster', num2str(id)];
        EasyPlot.exportFigure(fig, fullfile(output_folder, output_filename));
        close all;
    end

    % remove the outliers
    removeSpikes(folder_data, id, spike_ids(is_outliers));
    fprintf('Removed %d / %d spikes in cluster %d!\n', sum(is_outliers), length(spike_ids), id);

    if mod(k, 100) == 1
        fprintf('%d / %d done!\n', k, length(cluster_non_noise));
    end
end

end



























