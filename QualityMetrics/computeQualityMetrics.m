function computeQualityMetrics(folder_data, user_settings)
% COMPUTEQUALITYMETRICS compute all the quality metrics of non-noise clusters
%
% - Input
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
% 
% Output:
%   QualityMetrics.mat will be generated
%
% The default value from https://github.com/AllenInstitute/ecephys_spike_sorting/blob/919992748a5324724ba87169ecdcf9eb6e3b9973/ecephys_spike_sorting/scripts/create_input_json.py#L185
% "quality_metrics_params" : {
%     "isi_threshold" : 0.0015,
%     "min_isi" : 0.000166,
%     "num_channels_to_compare" : 7,
%     "max_spikes_for_unit" : 500,
%     "max_spikes_for_nn" : 10000,
%     "n_neighbors" : 4,
%     'n_silhouette' : 10000,
%     "quality_metrics_output_file" : os.path.join(kilosort_output_directory, "metrics_test.csv"),
%     "drift_metrics_interval_s" : 51,
%     "drift_metrics_min_spikes_per_interval" : 10,
% 
%     "include_pc_metrics" : True
% 
% From https://github.com/AllenInstitute/ecephys_spike_sorting/blob/919992748a5324724ba87169ecdcf9eb6e3b9973/ecephys_spike_sorting/modules/quality_metrics/_schemas.py#L7
%
% isi_threshold = Float(required=False, default=0.0015, help='Maximum time (in seconds) for ISI violation')
% min_isi = Float(required=False, default=0.00, help='Minimum time (in seconds) for ISI violation')
% num_channels_to_compare = Int(required=False, default=13, help='Number of channels to use for computing PC metrics; must be odd')
% max_spikes_for_unit = Int(required=False, default=500, help='Number of spikes to subsample for computing PC metrics')
% max_spikes_for_nn = Int(required=False, default=10000, help='Further subsampling for NearestNeighbor calculation')
% n_neighbors = Int(required=False, default=4, help='Number of neighbors to use for NearestNeighbor calculation')
% n_silhouette = Int(required=False, default=10000, help='Number of spikes to use for calculating silhouette score')
% 
% drift_metrics_min_spikes_per_interval = Int(required=False, default=10, help='Minimum number of spikes for computing depth')
% drift_metrics_interval_s = Float(required=False, default=100, help='Interval length is seconds for computing spike depth')
% 
% quality_metrics_output_file = String(required=True, help='CSV file where metrics will be saved')
% 
% include_pc_metrics = Bool(required=False, default=True, help='Compute features that require principal components')
%

if nargin < 1
    folder_data = './';
end

if nargin < 2
    num_channels_to_compare = 7;
    max_spikes_for_unit = 500;
    max_spikes_for_nn = 10000;
    n_neighbors = 4;
    n_silhouette = 10000;
else
    num_channels_to_compare = user_settings.num_channels_to_compare;
    max_spikes_for_unit = user_settings.max_spikes_for_unit;
    max_spikes_for_nn = user_settings.max_spikes_for_nn;
    n_neighbors = user_settings.n_neighbors;
    n_silhouette = user_settings.n_silhouette;    
end

% read the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));
cluster_ids = unique(spike_clusters);

n_cluster = length(cluster_ids);

cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
labels = cell(n_cluster, 1);
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    for k = 1:n_cluster
        idx = find(cluster_group.cluster_id == cluster_ids(k));
        if isempty(idx)
            continue
        end

        labels{k} = cluster_group.group{idx};
    end
end

% spike_locations = zeros(n_cluster, 2);
isi_violations = NaN(n_cluster, 1);
amplitude_cutoffs = NaN(n_cluster, 1);
presence_ratio = NaN(n_cluster, 1);
amplitude_median = NaN(n_cluster, 1);

disp('Computing Non-PC features...');
for k = 1:n_cluster
    if strcmpi(labels{k}, 'noise')
        if mod(k, 10) == 1
            fprintf('%d / %d done!\n', k, n_cluster);
        end
        continue
    end
    
    spike_time_this = double(spike_times(spike_clusters==cluster_ids(k)))./30000*1000; % in ms
    amplitude_this = double(amplitudes(spike_clusters==cluster_ids(k)));
    amplitude_this = rmoutliers(amplitude_this, 'median', 'ThresholdFactor', 5);

    t_begin = 0;
    t_end = double(max(spike_times))./30000*1000;

    isi_violations(k) = isiViolations(spike_time_this);
    amplitude_cutoffs(k) = amplitudeCutoffs(amplitude_this);
    amplitude_median(k) = median(amplitude_this);
    presence_ratio(k) = presenceRatio(spike_time_this, t_begin, t_end);
    
    if mod(k, 10) == 1
        fprintf('%d / %d done!\n', k, n_cluster);
    end
end

%% PC based metrics
% load files
chanMap = load(fullfile(folder_data, 'chanMap.mat'));
pc_features = readNPY(fullfile(folder_data, 'pc_features.npy'));
pc_feature_ind = double(readNPY(fullfile(folder_data, 'pc_feature_ind.npy')));
spike_templates = double(readNPY(fullfile(folder_data, 'spike_templates.npy')));

assert(mod(num_channels_to_compare, 2) == 1);
half_spread = floor((num_channels_to_compare - 1) / 2);

cluster_ids = double(unique(spike_clusters));
template_ids = double(unique(spike_templates));

template_peak_channels = zeros(length(template_ids), 1);
cluster_peak_channels = zeros(length(cluster_ids), 1); % 0 ~ n_channel-1

for idx = 1:length(template_ids)
    template_id = template_ids(idx);
    for_template = squeeze(spike_templates == template_id);
    pc_max = find(mean(pc_features(for_template, 1, :), 1) == max(mean(pc_features(for_template, 1, :), 1)), 1);
    template_peak_channels(idx) = pc_feature_ind(template_id+1, pc_max); % matlab index from 1
end

for idx = 1:length(cluster_ids)
    cluster_id = cluster_ids(idx);
    for_unit = squeeze(spike_clusters == cluster_id);
    templates_for_unit = unique(spike_templates(for_unit));
    template_positions = ismember(template_ids, templates_for_unit);
    cluster_peak_channels(idx) = round(median(template_peak_channels(template_positions)));
end

isolation_distance = NaN(n_cluster, 1);
d_prime = NaN(n_cluster, 1);
nn_miss_rate = NaN(n_cluster, 1);
nn_hit_rate = NaN(n_cluster, 1);
l_ratio = NaN(n_cluster, 1);

disp('Computing PC features...');
for idx = 1:n_cluster
    if strcmpi(labels{idx}, 'noise')
        continue
    end

    cluster_id = cluster_ids(idx);
    [isolation_distance(idx), d_prime(idx), nn_miss_rate(idx), nn_hit_rate(idx), l_ratio(idx)] = ...
        calculate_pc_metrics_one_cluster(...
        cluster_peak_channels, idx, cluster_id, cluster_ids,...
        half_spread, pc_features, pc_feature_ind,...
        spike_clusters, spike_templates,...
        max_spikes_for_unit, max_spikes_for_nn, n_neighbors, chanMap);

    if mod(idx, 10) == 1
        fprintf('%d / %d done!\n', idx, n_cluster);
    end
end

%% save to files
metric_names = {'ISI violations', 'Amplitude cutoffs', 'Presence ratio', 'Median Amplitude', 'Isolation distance', 'D prime', 'Nearest-neighbor miss rate', 'Nearest-neighbor hit rate', 'L ratio'};
metrics = {isi_violations, amplitude_cutoffs, presence_ratio, amplitude_median, isolation_distance, d_prime, nn_miss_rate, nn_hit_rate, l_ratio};

save(fullfile(folder_data, 'QualityMetrics.mat'),...
    'cluster_ids', 'isi_violations', 'amplitude_cutoffs', 'presence_ratio', 'labels', 'amplitude_median',...
    'isolation_distance', 'd_prime', 'nn_miss_rate', 'nn_hit_rate', 'l_ratio',...
    'metrics', 'metric_names');

end