function getPotentialSplits(folder_data, user_settings)
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
verbose = user_settings.merging.verbose;

% load the data
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));
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
split_info = {}; % cluster_id, spike_ids_A, spike_ids_B

for k = 1:length(cluster_id_non_noise)
    spike_ids = find(spike_clusters == cluster_id_non_noise(k));
    amplitudes_this = amplitudes(spike_ids);
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
    % do dbscan on the waveform
    step_size = 3;
    flag = false;
    for j = 1:n_channels_included
        for i = 1:2:size(waveforms, 3)
            data = squeeze(waveforms(:,ch_included(j),i));

            AIC = zeros(1,2);
            GMModels = cell(1,2);
            flag_fitted = true;
            for ii = 1:2
                try
                    GMModels{ii} = fitgmdist(data, ii, 'CovarianceType', 'diagonal');
                    AIC(ii)= GMModels{ii}.AIC;
                catch
                    flag_fitted = false;
                end
            end

            if ~flag_fitted
                continue
            end
            
            [~, numComponents] = min(AIC);
            overlap_ratio = overlapRatio(GMModels{2});
            JS_divergence = JSdivergence(GMModels{2}.mu(1), GMModels{2}.Sigma(1), GMModels{2}.mu(2), GMModels{2}.Sigma(2));

            if numComponents > 1 && GMModels{2}.Converged && overlap_ratio < 0.8
                overlap_ratio
                JS_divergence

                data_all = mmap.Data.x(ch_included(j), spike_times_this + waveform_window(1)+i-1);
                data_all = squeeze(double(data_all'));

                mdl = fitgmdist(data_all, numComponents, 'CovarianceType', 'diagonal');
                idx = mdl.cluster(data_all);
                if sum(idx==1) < length(spike_times_this)*0.1 || sum(idx==2) < length(spike_times_this)*0.1
                    continue
                else
                    flag = true;

                    figure;
                    histogram(data_all(idx==1));
                    hold on;
                    histogram(data_all(idx==2));
                    pause(1);
                    close all;

                    break
                end

               
            end
        end

        if flag
            break
        end
    end

    if ~flag
        fprintf('%d / %d done!\n', k, length(cluster_id_non_noise));
        continue
    end


    spike_ids_A = spike_ids(idx==1);
    spike_ids_B = spike_ids(idx==2);
    
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

end



