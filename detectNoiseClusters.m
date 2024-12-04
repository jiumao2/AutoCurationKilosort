function detectNoiseClusters(folder_data, user_settings)
% DETECTNOISECLUSTERS detect and remove the clusters which are pure noise
%
% Input:
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
%
% Output:
%   The dectected noise clusters will be labeled as "noise" in cluster_group.tsv
%
% Determination of noise:
%   (1) Firing rate < 0.05 Hz
%   (2) Signal-to-noise ratio > 2
%

% load the params
min_firing_rate = user_settings.detectNoiseClusters.min_firing_rate;
min_signal_to_noise_ratio = user_settings.detectNoiseClusters.min_signal_to_noise_ratio;
n_random_spikes = user_settings.detectNoiseClusters.n_random_spikes;
waveform_window = user_settings.detectNoiseClusters.waveform_window;
baseline_window = user_settings.detectNoiseClusters.baseline_window - waveform_window(1) + 1;

% load the data
path_data = fullfile(folder_data, 'temp_wh.dat');
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
load(fullfile(folder_data, 'ops.mat'));

dir_output = dir(path_data);
nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});

% detect noise
disp('Start detecting noise clusters');
duration_sec = double(max(spike_times) - min(spike_times))./30000;

cluster_ids = unique(spike_clusters);
idx_noise = [];
for k = 1:length(cluster_ids)
    id = cluster_ids(k);
    spike_ids = find(spike_clusters == id);

    fr = length(spike_ids)./duration_sec;
    if fr < min_firing_rate
        idx_noise = [idx_noise, id];
        continue
    end

    spike_times_this = spike_times(spike_ids);

    n_waveforms = min(length(spike_times_this), n_random_spikes);
    idx_rand = randperm(length(spike_times_this), n_waveforms);
    waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
    for j = 1:n_waveforms
        waveforms(j,:,:) = mmap.Data.x(:, spike_times_this(idx_rand(j)) + waveform_window(1):spike_times_this(idx_rand(j)) + waveform_window(2));
    end
    
    mean_waveforms = squeeze(mean(waveforms, 1)); % 383 x 64
    [amplitude, ch_largest] = max(max(mean_waveforms,[],2) - min(mean_waveforms,[],2));
    baselines = waveforms(:, ch_largest, baseline_window(1):baseline_window(2));
    variance_baseline = mean(baselines(:).^2);

    snr_this = amplitude.^2./variance_baseline;
    if snr_this < min_signal_to_noise_ratio
        fprintf('[snr = %.2f] Cluster %d is putative noise!\n', snr_this, id);
        idx_noise = [idx_noise, id];
        continue
    end
    
    if mod(k, 50) == 1
        fprintf('%d / %d done!\n', k, length(cluster_ids));
    end
end

% remove the noise clusters
fprintf('Found %d / %d noise clusters!\n', length(idx_noise), length(cluster_ids));
labelKilosort(folder_data, idx_noise, 'noise');

end
