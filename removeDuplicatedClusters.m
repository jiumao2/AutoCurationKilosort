function removeDuplicatedClusters(folder_data, user_settings)
% REMOVEDUPLICATEDCLUSTERS remove duplicated clusters which are actually from
% the same units but detected multiple times.
%
% Input:
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
%
% Output:
%   Removed clusters are labeled as "noise" in cluster_group.tsv
%
% Determination of duplicated clusters:
%   (1) Lots of very close spike times (dt <= 0.5, proportion > 10%)
%   (2) The same channel with the largest amplitude
% 
% The duplicated ones with less spike number will be removed. Noise
% clusters are excluded in this process.
%

n_random_spikes = user_settings.duplicatedClusters.n_random_spikes;
waveform_window = user_settings.duplicatedClusters.waveform_window;

% load the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
cluster_ids = unique(spike_clusters);

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
ch_largest = zeros(1, length(cluster_id_non_noise));

disp('Computing the channel with largest amplitudes...');
for k = 1:length(cluster_id_non_noise)
    spike_times_this = spike_times(spike_clusters == cluster_id_non_noise(k));
    n_waveforms = min(length(spike_times_this), n_random_spikes);
    idx_rand = randperm(length(spike_times_this), n_waveforms);
    waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64

    idx_remove = [];
    for j = 1:n_waveforms
        if spike_times_this(idx_rand(j)) + waveform_window(2) > size(mmap.Data.x, 2)
            idx_remove = [idx_remove, j];
            continue
        end
        waveforms(j,:,:) = mmap.Data.x(:,...
            spike_times_this(idx_rand(j)) + waveform_window(1):spike_times_this(idx_rand(j)) + waveform_window(2));
    end
    waveforms(idx_remove,:,:) = [];

    mean_waveforms = squeeze(mean(waveforms, 1)); % 383 x 64
    [~, ch_largest(k)] = max(max(mean_waveforms,[],2) - min(mean_waveforms,[],2));
    
    if mod(k, 50) == 1
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
    end
end

disp('Detecting duplicated clusters...');
for k = 1:length(cluster_id_non_noise)
    for j = k+1:length(cluster_id_non_noise)
        % Dulicated clusters should have same ch_largest
        if ch_largest(k) ~= ch_largest(j)
            continue
        end

        is_duplicated = isDuplicatedCluster(folder_data,...
            cluster_id_non_noise(k),...
            cluster_id_non_noise(j),...
            user_settings,...
            spike_times,...
            spike_clusters);

        if is_duplicated
            % remove the one with less spike number
            num_spikes_A = sum(spike_clusters == cluster_id_non_noise(k));
            num_spikes_B = sum(spike_clusters == cluster_id_non_noise(j));
            if num_spikes_A >= num_spikes_B
                labelKilosort(folder_data, cluster_id_non_noise(j), 'noise');
            else
                labelKilosort(folder_data, cluster_id_non_noise(k), 'noise');
            end
        end
    end
    
    if mod(k, 50) == 1
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
    end
end

end