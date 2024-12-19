folder_data = 'D:\Pierce\20241116\catgt_Exp_g0';
setting_filenames = 'C:\Users\jiumao\Desktop\AutoCurationKilosort\settings.json';

% read the settings
userSettings = jsonc.jsoncDecode(fileread(setting_filenames));

% remove clusters which are pure noise
detectNoiseClusters(folder_data, userSettings);

% clean the waveforms in each cluster
removeNoiseInsideCluster(folder_data, userSettings);

% todo: detect and do splits and merges
% split first (should be very conserative)
split_info = getPotentialSplits(folder_data, userSettings);

% merge
getPotentialMerges(folder_data, userSettings);

% remove duplicated clusters
removeDuplicatedClusters(folder_data, userSettings);

% determine the quality of each cluster
computeQualityMetrics(folder_data);
labelWithQualityMetrics(folder_data, userSettings);

% realign the spike times
realignClusterSpikeTimes(folder_data, userSettings);

% output to cluster_info.tsv
updateClusterInfo(folder_data);
