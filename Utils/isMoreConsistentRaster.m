function [out, info] = isMoreConsistentRaster(spike_times_A_ms, spike_times_B_ms, event_times_ms, user_settings)

% load the params
t_pre = user_settings.rasterMetric.t_pre_ms;
t_post = user_settings.rasterMetric.t_post_ms;
binwidth = user_settings.rasterMetric.binwidth_ms;
temporal_corr_thres = user_settings.rasterMetric.temporal_corr_thres;
peth_corr_thres = user_settings.rasterMetric.peth_corr_thres;

% compute the spike counts in each trial and time points
t_edges = t_pre:binwidth:t_post;
t_bins = 0.5*(t_edges(1:end-1)+t_edges(2:end));

n_trial = length(event_times_ms);
spike_times = {spike_times_A_ms, spike_times_B_ms};
spike_counts = {zeros(n_trial, length(t_bins)), zeros(n_trial, length(t_bins))};

for k = 1:n_trial
    t_event = event_times_ms(k);
    for j = 1:length(spike_times)
        st = spike_times{j} - t_event;
        for i = 1:length(t_bins)
            spike_counts{j}(k,i) = sum(st > t_edges(i) & st < t_edges(i+1));
        end
    end
end

% compute the mean of the spike counts (PETH)
mean_spike_counts = {mean(spike_counts{1}), mean(spike_counts{2})};
spike_number_per_trial = {mean(spike_counts{1},2), mean(spike_counts{2},2)};

% the spike number should be complementary
temp = corrcoef(spike_number_per_trial{1}, spike_number_per_trial{2});
temporal_corr = temp(1,2);
temp = corrcoef(mean_spike_counts{1}, mean_spike_counts{2});
peth_corr = temp(1,2);

info.temporal_corr = temporal_corr;
info.peth_corr = peth_corr;

if temporal_corr > temporal_corr_thres
    out = false;
    return
end

if peth_corr < peth_corr_thres
    out = false;
    return
end

out = true;

end