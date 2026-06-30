function n_channels = getKilosortNChannels(folder_data)
%GETKILOSORTNCHANNELS Return the effective Kilosort channel count.
%
% Prefer channel_map.npy because it is a numeric Kilosort output available
% in both Kilosort 2.5 and Kilosort 4. Fall back to ops.mat for older
% outputs that do not include channel_map.npy.

path_channel_map = fullfile(folder_data, 'channel_map.npy');
path_channel_positions = fullfile(folder_data, 'channel_positions.npy');
path_ops = fullfile(folder_data, 'ops.mat');

if isfile(path_channel_map)
    channel_map = readNPY(path_channel_map);
    n_channels = numel(channel_map);

    if isfile(path_channel_positions)
        channel_positions = readNPY(path_channel_positions);
        if size(channel_positions, 1) ~= n_channels
            error(['channel_positions.npy has %d rows, but channel_map.npy ', ...
                'contains %d channels in %s.'], ...
                size(channel_positions, 1), n_channels, folder_data);
        end
    end
    return
end

if isfile(path_ops)
    data = load(path_ops, 'ops');
    if isfield(data, 'ops') && isfield(data.ops, 'Nchan')
        n_channels = data.ops.Nchan;
        return
    end

    error('ops.mat exists but does not contain ops.Nchan in %s.', folder_data);
end

error(['Cannot determine Kilosort channel count in %s. Expected ', ...
    'channel_map.npy or ops.mat with ops.Nchan.'], folder_data);

end
