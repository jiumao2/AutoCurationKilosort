pairs = {[1799, 1871], [1715, 1717], [663, 1095], [518, 522]};
info = [];
for k = 1:length(pairs)
    for j = 1:4
        spike_times_A = spike_times(spike_clusters==pairs{k}(1));
        spike_times_B = spike_times(spike_clusters==pairs{k}(2));
    
        [~, info_this] = isMoreConsistentRaster(...
            double(spike_times_A)./30,...
            double(spike_times_B)./30,...
            event_times_ms{j},...
            user_settings);
    
        if isempty(info)
            info = info_this;
        else
            info(k,j) = info_this;
        end
    end
end

