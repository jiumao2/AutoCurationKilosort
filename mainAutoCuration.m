path_autocuration = 'path_to_AutoCurationKilosort\AutoCurationKilosort';
folder_data = './catgt_Exp_g0';
setting_filenames = 'path_to_settings/settings.json';

% add the path
addpath(path_autocuration);
addpath(fullfile(path_autocuration, 'PhyManipulations/'));
addpath(fullfile(path_autocuration, 'QualityMetrics/'));
addpath(genpath(fullfile(path_autocuration, 'Utils/')));

% read the settings
userSettings = jsonc.jsoncDecode(fileread(setting_filenames));

% remove clusters which are pure noise
detectNoiseClusters(folder_data, userSettings);

% clean the waveforms in each cluster
removeNoiseInsideCluster(folder_data, userSettings);

% remove duplicated clusters
removeDuplicatedClusters(folder_data, userSettings);

% determine the quality of each cluster
computeQualityMetrics(folder_data, userSettings);
labelWithQualityMetrics(folder_data, userSettings);

% realign the spike times
realignClusterSpikeTimes(folder_data, userSettings);

% output to cluster_info.tsv
updateClusterInfo(folder_data);
