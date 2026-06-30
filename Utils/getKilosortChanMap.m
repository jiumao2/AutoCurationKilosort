function chanMap = getKilosortChanMap(folder_data)
%GETKILOSORTCHANMAP Return channel geometry in chanMap.mat-compatible form.

path_chan_map_mat = fullfile(folder_data, 'chanMap.mat');
path_channel_map = fullfile(folder_data, 'channel_map.npy');
path_channel_positions = fullfile(folder_data, 'channel_positions.npy');

if isfile(path_chan_map_mat)
    chanMap = load(path_chan_map_mat);
    return
end

if ~isfile(path_channel_map) || ~isfile(path_channel_positions)
    error(['Cannot determine channel geometry in %s. Expected chanMap.mat ', ...
        'or both channel_map.npy and channel_positions.npy.'], folder_data);
end

channel_map = readNPY(path_channel_map);
channel_positions = readNPY(path_channel_positions);

if size(channel_positions, 1) ~= numel(channel_map)
    error(['channel_positions.npy has %d rows, but channel_map.npy ', ...
        'contains %d channels in %s.'], ...
        size(channel_positions, 1), numel(channel_map), folder_data);
end

chanMap.chanMap = double(channel_map(:)) + 1;
chanMap.connected = true(numel(channel_map), 1);
chanMap.xcoords = double(channel_positions(:, 1));
chanMap.ycoords = double(channel_positions(:, 2));

end
