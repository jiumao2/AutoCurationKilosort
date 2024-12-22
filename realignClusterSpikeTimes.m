function realignClusterSpikeTimes(folder_data, user_settings)
% REALIGNSPIKETIMES realign the spike times of good/mua clusters to center
% the peaks / troughs. It is helpful for getting a better mean waveforms
% and tracking units across days.
%
% Input:
%   - folder_data: the folder where the data is located
%   - user_settings: the global settings
%
% Output:
%   The modified spike times will be saved in spike_times.npy
%

n_random_spikes = user_settings.realignSpikeTimes.n_random_spikes;
waveform_window = user_settings.realignSpikeTimes.waveform_window;
baseline_window = user_settings.realignSpikeTimes.baseline_window - waveform_window(1) + 1;
n_channels_included = user_settings.realignSpikeTimes.n_channels_included;
verbose = user_settings.realignSpikeTimes.verbose;

% load the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');

path_data = fullfile(folder_data, 'temp_wh.dat');
load(fullfile(folder_data, 'ops.mat'));

dir_output = dir(path_data);
nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});

% get the good and mua clusters
cluster_ids = [];
labels = {};
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    cluster_ids = cluster_group.cluster_id(strcmpi(cluster_group.group, 'good') | strcmpi(cluster_group.group, 'mua'));
    labels = cluster_group.group(strcmpi(cluster_group.group, 'good') | strcmpi(cluster_group.group, 'mua'));
end

% realign
disp('Start realigning the spikes!');
idx_center_raw = -waveform_window(1)+1;

for k = 1:length(cluster_ids)
    spike_ids = find(spike_clusters == cluster_ids(k));
    spike_times_this = spike_times(spike_ids);

    % get the largest channel
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
    [~, idx_sort] = sort(max(mean_waveforms,[],2) - min(mean_waveforms,[],2), 'descend');
    ch_largest = idx_sort(1);
    
    baseline = mean(mean_waveforms(baseline_window));
    if abs(max(mean_waveforms(ch_largest,:))- baseline) < abs(min(mean_waveforms(ch_largest,:))- baseline)
        unit_type = 1;
        [~, idx_min] = min(mean_waveforms(ch_largest,:));
        dsample_template = idx_min - idx_center_raw;
    else
        unit_type = 2;
        [~, idx_max] = max(mean_waveforms(ch_largest,:));
        dsample_template = idx_max - idx_center_raw;
    end

    % realign the spike times
    assert(length(dsample_template) == 1);
    spike_times_this = spike_times_this + uint64(dsample_template);
    spike_times(spike_ids) = spike_times_this;
    
    if mod(k, 20) == 1
        fprintf('%d / %d done!\n', k, length(cluster_ids));
    end

    if verbose
        fig = EasyPlot.figure();
        ax = EasyPlot.axes(fig,...
            'Width', 3,...
            'Height', 3,...
            'MarginLeft', 0.1,...
            'MarginRight', 0.1,...
            'XAxisVisible','on',...
            'YAxisVisible', 'off');

        n_plot = min(100, n_waveforms);
        idx_plot = randperm(size(waveforms, 1), n_plot);

        x_waveform = waveform_window(1):waveform_window(2);
        plot(ax, x_waveform, squeeze(waveforms(idx_plot,ch_largest,:)), 'b-');

       
        x_plot = [];
        y_plot = [];
        for j = 1:length(idx_plot)
            x_plot = [x_plot, x_waveform-dsample_template, NaN];
            y_plot = [y_plot, squeeze(waveforms(idx_plot(j),ch_largest,:))', NaN];
        end
        plot(ax, x_plot, y_plot, 'r-');

        EasyPlot.setXLim(ax, [waveform_window(1), waveform_window(2)]);
        EasyPlot.setGeneralTitle(ax, ['Cluster ', num2str(cluster_ids(k)), ' (', labels{k}, ')']);
        EasyPlot.cropFigure(fig);

        output_folder = fullfile(folder_data, 'Fig/RealignClusterTemplates/');
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        output_filename = ['Cluster', num2str(cluster_ids(k))];
        EasyPlot.exportFigure(fig, fullfile(output_folder, output_filename));
        close all;
    end
end

% save the realigned spike times
updateSpikeTimes(folder_data, spike_times);
disp('Spike times are realigned!');

end