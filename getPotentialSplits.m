function split_info = getPotentialSplits(folder_data, user_settings)
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
min_cluster_percentage = user_settings.splitting.min_cluster_percentage;
n_random_spikes = user_settings.splitting.n_random_spikes;
waveform_window = user_settings.splitting.waveform_window;
n_channels_included = user_settings.splitting.n_channels_included;
hartigans_dip_test_alpha = user_settings.splitting.hartigans_dip_test_alpha;
verbose = user_settings.merging.verbose;

% load the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
path_data = fullfile(folder_data, 'temp_wh.dat');
load(fullfile(folder_data, 'ops.mat'));
dir_output = dir(path_data);
nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});

cluster_ids = unique(spike_clusters);

% get the non-noise clusters
cluster_id_noise = [];
cluster_group = readtable(fullfile(folder_data, 'cluster_group.tsv'), 'Delimiter', '\t', 'FileType', 'text');
if any(strcmpi(cluster_group.Properties.VariableNames, 'group'))
    cluster_id_noise = cluster_group.cluster_id(strcmpi(cluster_group.group, 'noise'));
end

cluster_id_non_noise = setdiff(cluster_ids, cluster_id_noise);


%% dectect potential splits with dbscan
disp('Start detecting potential splits!');
split_info = {}; % cluster_id, spike_ids_A, spike_ids_B

for k = 1:length(cluster_id_non_noise)
    spike_ids = find(spike_clusters == cluster_id_non_noise(k));
    spike_times_this = spike_times(spike_ids);

    n_waveforms = min(length(spike_times_this), n_random_spikes);
    idx_rand = randperm(length(spike_times_this), n_waveforms);
    waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
    for j = 1:n_waveforms
        waveforms(j,:,:) = mmap.Data.x(:,...
            spike_times_this(idx_rand(j)) + waveform_window(1):spike_times_this(idx_rand(j)) + waveform_window(2));
    end

    mean_waveforms = squeeze(mean(waveforms, 1)); % 383 x 64
    [~, idx_sort] = sort(max(mean_waveforms,[],2) - min(mean_waveforms,[],2), 'descend');
    
    ch_included = idx_sort(1:n_channels_included);

    % do HartigansDipSignifTest on each time point on the waveform
    dip_max = -1;
    split_place = [];
    for j = 1:n_channels_included
        for i = 1:size(waveforms, 3)
            data = squeeze(waveforms(:,ch_included(j),i));
            [h, dip] = HartigansDipSignifTest(data, max(0.1, hartigans_dip_test_alpha)); 

            if h==0 && dip > dip_max
                dip_max = dip;
                split_place = [j,i];
            end
        end
    end
    
    % continue if unimodal on every time point
    if dip_max == -1
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
        continue
    end
    
    % read the full data for splitting
    data_all = mmap.Data.x(ch_included(split_place(1)),...
        spike_times_this + waveform_window(1)+split_place(2)-1);
    data_all = squeeze(double(data_all'));

    idx = kmeans(data_all, 2);
    h = HartigansDipSignifTest(data_all, hartigans_dip_test_alpha);

    if h==1 || min(sum(idx==1), sum(idx==2)) < length(spike_times_this)*min_cluster_percentage/100
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
        continue
    end

    spike_ids_A = spike_ids(idx==1);
    spike_ids_B = spike_ids(idx==2);

    % to do: the splits should not pass the merging criteria (control the false positive rate)
    
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

fprintf('Found %d potential splits!\n', length(split_info));

end



