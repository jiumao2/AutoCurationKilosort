function is_duplicated = isDuplicatedCluster(folder_data, cluster_id_A, cluster_id_B, user_settings, spike_times, spike_clusters)
% ISDUPLICATEDCLUSTER check whether the two clusters are identical
%
% Input:
%   - folder_data: the folder where the data is located
%   - cluster_id_A: the first cluster ID
%   - cluster_id_B: the other cluster ID
%   - user_settings: the global settings
%   - spike_times: optional. It would be faster if given rather than
%   reading from data.
%   - spike_clusters: optional. It would be faster if given rather than
%   reading from data.
%
% Output:
%   - is_duplicated: true or false. 
%
% Determination of duplicated clusters:
%   (1) Lots of very close spike times (dt <= 0.5, proportion > 10%)
%   (2) The same channel with the largest amplitude
%

% load the params
dt_sample = user_settings.duplicatedClusters.dt * 30;
overlap_percentage = user_settings.duplicatedClusters.overlap_percentage;
n_random_spikes = user_settings.duplicatedClusters.n_random_spikes;
waveform_window = user_settings.duplicatedClusters.waveform_window;
verbose = user_settings.duplicatedClusters.verbose;

% load the data if not given
if nargin < 5
    spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
end
if nargin < 6
    spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
end

stA = int64(spike_times(spike_clusters == cluster_id_A));
stB = int64(spike_times(spike_clusters == cluster_id_B));

if isempty(stA) || isempty(stB)
    is_duplicated = false;
    return
end

% compute the percentage of overlap
idx_nearest = findNearestPoint(stB, stA);
num_overlap_A = sum(abs(stA - stB(idx_nearest)) < dt_sample);

idx_nearest = findNearestPoint(stA, stB);
num_overlap_B = sum(abs(stB - stA(idx_nearest)) < dt_sample);

if num_overlap_A / length(stA) < overlap_percentage/100 && num_overlap_B / length(stB) < overlap_percentage/100
    is_duplicated = false;
    return
end

% check the channel with the largest amplitude
path_data = fullfile(folder_data, 'temp_wh.dat');
load(fullfile(folder_data, 'ops.mat'));
dir_output = dir(path_data);
nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});

ch_largest = zeros(1, 2);
st_all = {stA, stB};
waveforms_all = cell(1,2);

for k = 1:2
    n_waveforms = min(length(st_all{k}), n_random_spikes);
    idx_rand = randperm(length(st_all{k}), n_waveforms);
    waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
    for j = 1:n_waveforms
        waveforms(j,:,:) = mmap.Data.x(:,...
            st_all{k}(idx_rand(j)) + waveform_window(1):st_all{k}(idx_rand(j)) + waveform_window(2));
    end

    waveforms_all{k} = waveforms;

    mean_waveforms = squeeze(mean(waveforms, 1)); % 383 x 64
    [~, ch_largest(k)] = max(max(mean_waveforms,[],2) - min(mean_waveforms,[],2));
end

if ch_largest(1) == ch_largest(2)
    is_duplicated = true;
    outcome_str = 'duplicated';
else
    is_duplicated = false;
    outcome_str = 'different';
end

fprintf('Cluster %d  and cluster %d are %s!\n', cluster_id_A, cluster_id_B, outcome_str);

% make the figure for visualization
if verbose
    fig = EasyPlot.figure();
    
    ax_corr = EasyPlot.axes(fig,...
        'Width', 3,...
        'Height', 3,...
        'YAxisVisible', 'off',...
        'MarginBottom', 1,...
        'MarginLeft', 0.1,...
        'MarginRight', 0.1);
    ax_waveforms_A = EasyPlot.createAxesAgainstAxes(fig, ax_corr, 'right',...
        'XAxisVisible', 'off',...
        'YAxisVisible', 'off');
    ax_waveforms_B = EasyPlot.createAxesAgainstAxes(fig, ax_waveforms_A, 'right',...
        'XAxisVisible', 'off',...
        'YAxisVisible', 'off');
    
    % cross correlogram
    s = cell(1,2);
    max_t = max(max(stA ./ 30), max(stB ./ 30));
    for k = 1:length(s)
        temp = zeros(1, max_t);
        temp(st_all{k} ./ 30) = 1;
        s{k} = temp;
    end

    [r, lags] = xcorr(s{1}, s{2}, 50);
    r(lags == 0) = 0;
    bar(ax_corr, lags, r, 'black');
    xlim(ax_corr, [-50, 50]);

    % waveforms in the ch_largest
    plot(ax_waveforms_A,...
        waveform_window(1):waveform_window(2),...
        squeeze(waveforms_all{1}(:, ch_largest(1), :)), 'k-');
    plot(ax_waveforms_A,...
        waveform_window(1):waveform_window(2),...
        squeeze(waveforms_all{2}(:, ch_largest(1), :)), 'r-');

    p1 = plot(ax_waveforms_B,...
        waveform_window(1):waveform_window(2),...
        squeeze(waveforms_all{2}(:, ch_largest(2), :)), 'r-');
    p2 = plot(ax_waveforms_B,...
        waveform_window(1):waveform_window(2),...
        squeeze(waveforms_all{1}(:, ch_largest(2), :)), 'k-');

    EasyPlot.legend(ax_waveforms_B,...
        {['Cluster ', num2str(cluster_id_B)], ['Cluster ', num2str(cluster_id_A)]},...
        'selectedPlots', [p1(1), p2(1)],...
        'location', 'northeastoutside');

    EasyPlot.setXLim({ax_waveforms_A, ax_waveforms_B}, waveform_window);
    EasyPlot.setYLim({ax_waveforms_A, ax_waveforms_B});
    title(ax_waveforms_A, ['Ch ', num2str(ch_largest(1))]);
    title(ax_waveforms_B, ['Ch ', num2str(ch_largest(2))]);

    if is_duplicated
        outcome_str = 'duplicated';
    else
        outcome_str = 'different';
    end

    h = EasyPlot.setGeneralTitle({ax_corr, ax_waveforms_A, ax_waveforms_B},...
        ['Cluster#' num2str(cluster_id_A), ' (n = ', num2str(length(stA)), ')',...
        ' VS Cluster#', num2str(cluster_id_B), ' (n = ', num2str(length(stB)), ')',...
        ', ', outcome_str]);
    EasyPlot.move(h, 'dy', 0.3);

    output_folder = fullfile(folder_data, 'Fig/DuplicatedClusters/');
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    output_filename = ['Cluster', num2str(cluster_id_A), '_VS_', 'Cluster', num2str(cluster_id_B), '_', outcome_str];
    
    EasyPlot.cropFigure(fig);
    EasyPlot.exportFigure(fig, fullfile(output_folder, output_filename), 'dpi', 300);
    close all;
end


end