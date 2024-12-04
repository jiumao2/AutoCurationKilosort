function updateSpikeTimes(folder_data, spike_times_new)
% UPDATESPIKETIMES update the spike times and the corresponding files
%
% Input:
%   - folder_data: 
%   - spike_times_new
%
% Output:
%   the files to be update:
%       (1) spike_times.npy
%       (2) spike_clusters.npy
%       (3) amplitudes.npy
%       (4) pc_features.npy
%       (5) spike_templates.npy
%       (6) template_features.npy
%

% read all npy files
amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
spike_templates = readNPY(fullfile(folder_data, 'spike_templates.npy'));
pc_features = readNPY(fullfile(folder_data, 'pc_features.npy'));
template_features = readNPY(fullfile(folder_data, 'template_features.npy'));

% the spike times should be increasing values
assert(length(spike_times) == length(spike_times_new));
[~, idx_sort] = sort(spike_times_new);

writeNPY(spike_times_new(idx_sort), fullfile(folder_data, 'spike_times.npy'));
writeNPY(amplitudes(idx_sort), fullfile(folder_data, 'amplitudes.npy'));
writeNPY(spike_clusters(idx_sort), fullfile(folder_data, 'spike_clusters.npy'));
writeNPY(spike_templates(idx_sort), fullfile(folder_data, 'spike_templates.npy'));
writeNPY(pc_features(idx_sort, :, :), fullfile(folder_data, 'pc_features.npy'));
writeNPY(template_features(idx_sort, :), fullfile(folder_data, 'template_features.npy'));
end