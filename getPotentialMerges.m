function merge_pairs = getPotentialMerges(folder_data, user_settings)
% GETPOTENTIALMERGES find potential merges
%
% Input:
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
%
% Output:
%   To do
%
% Determination of potential merges:
%   (1) Exclude all noise clusters
%   (2) Do not merge duplicated clusters
%   (3) Do not merges clusters that are far away
%   (4) Merging increases the consistency of firing rate (variance of firing rate)
%   (5) Merging increases the consistency of PETH:
%           similar in PETH correlation and anticorrelated in the total spike
%           number in each trial
%   (6) Merging does not increase ISI violation very much
%   (7) If overlapped, the cross correlogram should be significantly different than a flat
%   line
%   (8) If overlapped, the -1~1 ms bin should be low in the cross
%   correlogram (sum of these bins should be smaller than mean cross correlation)
%

% read the params
n_random_spikes = user_settings.merging.n_random_spikes;
waveform_window = user_settings.merging.waveform_window;
max_distance_um = user_settings.merging.max_distance_um;
binwidth_firing_rate_sec = user_settings.merging.binwidth_firing_rate_sec;
% isi_violation_thres = user_settings.merging.isi_violation_thres;
min_mean_cross_corr = user_settings.merging.min_mean_cross_corr;
cross_corr_range = user_settings.merging.cross_corr_range;
cross_corr_gaussian_kernel = user_settings.merging.cross_corr_gaussian_kernel;
isRasterMetric = user_settings.merging.rasterMetric;
eventsIncluded = user_settings.merging.eventsIncluded;

verbose = user_settings.merging.verbose;

% load the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
channel_positions = readNPY(fullfile(folder_data, 'channel_positions.npy'));
channel_map_ind0 = readNPY(fullfile(folder_data, 'channel_map.npy'));

cluster_ids = unique(spike_clusters);
min_t_sec = double(min(spike_times)) ./ 30000;
max_t_sec = double(max(spike_times)) ./ 30000;

% get the non-noise clusters
cluster_id_noise = [];
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    cluster_id_noise = cluster_group.cluster_id(strcmpi(cluster_group.group, 'noise'));
end

cluster_id_non_noise = setdiff(cluster_ids, cluster_id_noise);

% get the channels of the largest amplitudes of each clusters
path_data = fullfile(folder_data, 'temp_wh.dat');
load(fullfile(folder_data, 'ops.mat'));
dir_output = dir(path_data);
nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});
depth = zeros(1, length(cluster_id_non_noise));

disp('Computing the position of units with the channel with largest amplitudes...');
for k = 1:length(cluster_id_non_noise)
    spike_times_this = spike_times(spike_clusters == cluster_id_non_noise(k));
    n_waveforms = min(length(spike_times_this), n_random_spikes);
    idx_rand = randperm(length(spike_times_this), n_waveforms);
    waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
    for j = 1:n_waveforms
        waveforms(j,:,:) = mmap.Data.x(:,...
            spike_times_this(idx_rand(j)) + waveform_window(1):spike_times_this(idx_rand(j)) + waveform_window(2));
    end

    mean_waveforms = squeeze(mean(waveforms, 1)); % 383 x 64
    [~, ch_largest] = max(max(mean_waveforms,[],2) - min(mean_waveforms,[],2));
    depth(k) = channel_positions(channel_map_ind0 == ch_largest-1, 2);
    
    if mod(k, 50) == 1
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
    end
end

% load event times
if isRasterMetric
    event_labels = readcell(fullfile(folder_data, 'event_labels.csv'));
    event_data = readcell(fullfile(folder_data, 'events.csv'));
    event_times_ms = cell(1, length(event_labels));
    n_event = length(event_labels);

    for i = 1:n_event
        for ii = 1:size(event_data, 2)
            data_this = event_data{i,ii};
            if ~ismissing(data_this)
                event_times_ms{i} = [event_times_ms{i}, data_this*1000];
            end
        end
    end
end

%% detect the potential merges
disp('Start detecting potential merges!');
merge_pairs = [];
for k = 1:length(cluster_id_non_noise)
    for j = k+1:length(cluster_id_non_noise)
        % the position of units should be close enough
        if abs(depth(k) - depth(j)) > max_distance_um
            if mod(k, 10) == 1
                fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
            end
            continue
        end

        % do not merge duplicated clusters
        if isDuplicatedCluster(folder_data,...
                cluster_id_non_noise(k),...
                cluster_id_non_noise(j),...
                user_settings,...
                spike_times,...
                spike_clusters)
            continue
        end

        % consistency of firing rate
        t_edges_sec = min_t_sec:binwidth_firing_rate_sec:max_t_sec;

        spike_times_A = spike_times(spike_clusters==cluster_id_non_noise(k));
        spike_times_B = spike_times(spike_clusters==cluster_id_non_noise(j));

        spike_train_A = zeros(1, length(t_edges_sec)-1);
        spike_train_B = zeros(1, length(t_edges_sec)-1);
        spike_train_merged = zeros(1, length(t_edges_sec)-1);
        for i = 1:length(t_edges_sec)-1
            t_start_sample = t_edges_sec(i)*30000;
            t_end_sample = t_edges_sec(i+1)*30000;
            spike_train_A(i) = sum(spike_times_A > t_start_sample & spike_times_A < t_end_sample);
            spike_train_B(i) = sum(spike_times_B > t_start_sample & spike_times_B < t_end_sample);
            spike_train_merged(i) = spike_train_A(i) + spike_train_B(i);
        end

        if var(spike_train_merged) >...
                min([var(spike_train_A), var(spike_train_B)])
            fprintf('[Cluster %d | Cluster %d] Failed to pass the timing complementarity test!\n',...
                    cluster_id_non_noise(k),...
                    cluster_id_non_noise(j));
            if mod(k, 10) == 1
                fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
            end
            continue
        end

        % check cross correlogram
        spike_times_all = {spike_times(spike_clusters==cluster_id_non_noise(k)),...
            spike_times(spike_clusters==cluster_id_non_noise(j))};
        s = cell(1,2);

        min_t = min(min(spike_times_all{1} ./ 30), min(spike_times_all{2} ./ 30));
        max_t = max(max(spike_times_all{1} ./ 30), max(spike_times_all{2} ./ 30));
        t_bins = min_t:max_t;
        for i = 1:length(s)
            temp = zeros(1, length(t_bins));
            temp(spike_times_all{i}./30 -min_t + 1) = 1;
            s{i} = temp;
        end

        [r, lags] = xcorr(s{1}, s{2}, cross_corr_range(2));

        r_this1 = r(lags >= cross_corr_range(1) & lags <= cross_corr_range(2));
        r_smoothed1 = smoothdata(r_this1, 'gaussian', cross_corr_gaussian_kernel*5);
        
        r_this2 = r(lags <= -cross_corr_range(1) & lags >= -cross_corr_range(2));
        r_smoothed2 = smoothdata(r_this2, 'gaussian', cross_corr_gaussian_kernel*5);
        
        r_this = [r_this1, r_this2];
        r_smoothed = [r_smoothed1, r_smoothed2];

        r_mean = mean(r_this);
        
        if r_mean > min_mean_cross_corr
            if any(r(lags>=-2 & lags<=2) > mean(r))
                fprintf('[Cluster %d | Cluster %d] Failed to pass the refactory period test!\n',...
                    cluster_id_non_noise(k),...
                    cluster_id_non_noise(j));
                if mod(k, 10) == 1
                    fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
                end
                continue
            end

            chi2_stat = sum((r_smoothed - r_mean).^2./r_mean);
            df = length(r_this) - 1;
            chi2_crit = chi2inv(1 - 0.05, df);
            
            if chi2_stat < chi2_crit
                fprintf('[Cluster %d | Cluster %d] Failed to pass cross correlogram metric!\n',...
                    cluster_id_non_noise(k),...
                    cluster_id_non_noise(j));
                if mod(k, 10) == 1
                    fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
                end
                continue
            end
        end

        % check the raster metrics
        if isRasterMetric
            is_improved = false(1, length(eventsIncluded));
            for i = 1:length(eventsIncluded)
                idx = strcmpi(event_labels, eventsIncluded{i});
                is_improved(i) = isMoreConsistentRaster(...
                    double(spike_times_A)./30,...
                    double(spike_times_B)./30,...
                    event_times_ms{idx},...
                    user_settings);
            end

            if ~any(is_improved)
                fprintf('[Cluster %d | Cluster %d] Failed to pass raster metric!\n',...
                    cluster_id_non_noise(k),...
                    cluster_id_non_noise(j));
                if mod(k, 10) == 1
                    fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
                end
                continue
            else
                fprintf('[Cluster %d | Cluster %d] Passed the final raster metric [%s]!\n',...
                    cluster_id_non_noise(k),...
                    cluster_id_non_noise(j),...
                    strjoin(eventsIncluded(is_improved), ', '));
            end
        end

        merge_pairs = [merge_pairs; cluster_id_non_noise(k), cluster_id_non_noise(j)];
        
        if verbose
            title_name = ['Cluster#' num2str(cluster_id_non_noise(k)),...
                ' (n = ', num2str(length(spike_times_A)), ')',...
                ' VS Cluster#', num2str(cluster_id_non_noise(j)),...
                ' (n = ', num2str(length(spike_times_B)), ')'];
            save_filename = fullfile(folder_data, 'Fig/PotentialMerges/',...
                ['Cluster' num2str(cluster_id_non_noise(k)),...
                '_VS_Cluster', num2str(cluster_id_non_noise(j))]);
            
            plotClusterComparison(folder_data,...
                find(spike_clusters==cluster_id_non_noise(k)),...
                find(spike_clusters==cluster_id_non_noise(j)),...
                title_name,...
                save_filename);
        end
    
    end
    if mod(k, 10) == 1
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
    end
end

fprintf('Found %d potential merges!\n', size(merge_pairs, 1));

end