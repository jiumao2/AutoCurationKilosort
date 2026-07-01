function n_channels = getKilosortNChannels(folder_data)
%GETKILOSORTNCHANNELS Return the binary channel count for temp_wh.dat.
%
% Kilosort/Phy store the binary row count in params.py as n_channels_dat.
% This can be larger than channel_map.npy when the binary keeps an extra
% sync channel. Fall back to ops.mat for older outputs.

path_params = fullfile(folder_data, 'params.py');
path_ops = fullfile(folder_data, 'ops.mat');
path_temp_wh = fullfile(folder_data, 'temp_wh.dat');

n_channels = [];

if isfile(path_params)
    params_text = fileread(path_params);
    token = regexp(params_text, 'n_channels_dat\s*=\s*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        n_channels = str2double(token{1});
    end
end

if isempty(n_channels) && isfile(path_ops)
    data = load(path_ops, 'ops');
    if isfield(data, 'ops') && isfield(data.ops, 'Nchan')
        n_channels = data.ops.Nchan;
    else
        error('ops.mat exists but does not contain ops.Nchan in %s.', folder_data);
    end
end

if isempty(n_channels)
    error(['Cannot determine Kilosort channel count in %s. Expected ', ...
        'params.py with n_channels_dat or ops.mat with ops.Nchan.'], ...
        folder_data);
end

if ~isscalar(n_channels) || ~isfinite(n_channels) || n_channels < 1 || ...
        n_channels ~= round(n_channels)
    error('Invalid Kilosort channel count %g in %s.', n_channels, folder_data);
end

if isfile(path_temp_wh)
    file_info = dir(path_temp_wh);
    if mod(file_info.bytes, 2 * n_channels) ~= 0
        error(['temp_wh.dat size in %s is not divisible by int16 x ', ...
            '%d channels.'], folder_data, n_channels);
    end
end

end
