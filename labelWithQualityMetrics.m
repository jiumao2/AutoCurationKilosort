function labelWithQualityMetrics(folder_data, user_settings)
% LABELWITHQUALITYMETRICS Labels each cluster with the criteria in
% settings.json and the computed QualityMetrics.mat
%
% Input:
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
%
% Output:
%   Labelled clusters which can be found in cluster_group.tsv
%

params = user_settings.qualityMetrics;

% Check the existance of QualityMetrics.mat
filename_QM = fullfile(folder_data, 'QualityMetrics.mat');
if ~exist(filename_QM, 'file')
    error([filename_QM, ' not found! Please run "computeQualityMetrics" first']);
end

QM = load(filename_QM);

% Plot the metrics
plotQualityMetrics(QM, fullfile(folder_data, 'Fig'));

% Check the cluster_ids in the quality metrics are the same as the data
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
cluster_ids_data = unique(spike_clusters);

assert(length(cluster_ids_data) == length(QM.cluster_ids));
assert(length(intersect(cluster_ids_data, QM.cluster_ids)) == length(cluster_ids_data));

% get good clusters
metric_fieldnames = fieldnames(params.goodCriteria);
metric_thres = zeros(1, length(metric_fieldnames));
for k = 1:length(metric_thres)
    metric_thres(k) = params.goodCriteria.(metric_fieldnames{k});
end

idx_name = cellfun(@(name)find(strcmpi(params.metricFieldNames, name)), metric_fieldnames);
metric_names = params.metricNames(idx_name);
metric_types = params.metricTypes(idx_name);
idx_good = [];

disp('Labeling the clusters using quality metrics!');
for k = 1:length(QM.cluster_ids)
    is_passed = true;
    for j = 1:length(metric_names)
        idx_metric = find(strcmpi(QM.metric_names, metric_names{j}));
        if isnan(QM.metrics{idx_metric}(k))
            is_passed = false;
            break
        end

        if strcmpi(metric_types{j}, 'greater') && QM.metrics{idx_metric}(k) < metric_thres(j)
            is_passed = false;
            break
        elseif strcmpi(metric_types{j}, 'less') && QM.metrics{idx_metric}(k) > metric_thres(j)
            is_passed = false;
            break
        end
    end
    
    if is_passed
        idx_good = [idx_good, QM.cluster_ids(k)];
    end
end

% get mua clusters
metric_fieldnames = fieldnames(params.muaCriteria);
metric_thres = zeros(1, length(metric_fieldnames));
for k = 1:length(metric_thres)
    metric_thres(k) = params.muaCriteria.(metric_fieldnames{k});
end

idx_name = cellfun(@(name)find(strcmpi(params.metricFieldNames, name)), metric_fieldnames);
metric_names = params.metricNames(idx_name);
metric_types = params.metricTypes(idx_name);
idx_mua = [];

for k = 1:length(QM.cluster_ids)
    is_passed = true;
    for j = 1:length(metric_names)
        idx_metric = find(strcmpi(QM.metric_names, metric_names{j}));
        if isnan(QM.metrics{idx_metric}(k))
            is_passed = false;
            break
        end

        if strcmpi(metric_types{j}, 'greater') && QM.metrics{idx_metric}(k) < metric_thres(j)
            is_passed = false;
            break
        elseif strcmpi(metric_types{j}, 'less') && QM.metrics{idx_metric}(k) > metric_thres(j)
            is_passed = false;
            break
        end
    end
    
    if is_passed
        idx_mua = [idx_mua, QM.cluster_ids(k)];
    end
end

idx_mua = setdiff(idx_mua, idx_good);

% Do not contain any units that has been labeled as "noise"
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
idx_noise = [];
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    idx_noise = cluster_group.cluster_id(strcmpi(cluster_group.group, 'noise'));
end

idx_good = setdiff(idx_good, idx_noise);
idx_mua = setdiff(idx_mua, idx_noise);

% check idx_good and idx_mua
assert(isempty(intersect(idx_good, idx_mua)));
fprintf('Found %d good units and %d multi units!\n', length(idx_good), length(idx_mua));

% % label these clusters
labelKilosort(folder_data, idx_good, 'good');
labelKilosort(folder_data, idx_mua, 'mua');

end