min_cluster_percentage = user_settings.splitting.min_cluster_percentage;
spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));
%%
cluster_id = 1913;

spike_ids = find(spike_clusters == cluster_id);
amplitudes_this = amplitudes(spike_ids);
spike_times_sec = double(spike_times(spike_ids))./30000;

min_num_spikes = min_cluster_percentage*length(spike_times_sec)/100;

mdl = TimeVaryingGaussian();
mdl.fit(amplitudes_this, spike_times_sec, 60*5);
p = mdl.pdf(amplitudes_this, spike_times_sec);
p_crit = normpdf(4);
idx = p < p_crit;

figure;
plot(spike_times_sec(idx==0), amplitudes_this(idx==0), 'b.');
hold on;
plot(spike_times_sec(idx==1), amplitudes_this(idx==1), 'r.');
plot(mdl.t_bins, mdl.mu, 'k-', 'LineWidth', 2);




