function path_data = getKilosortTempWhPath(folder_data)
%GETKILOSORTTEMPWHPATH Return the preprocessed Kilosort waveform file path.
%
% The current waveform pipeline reads temp_wh.dat. It does not silently
% switch to raw .ap.bin because raw files may include sync channels and have
% different filtering/whitening semantics.

path_data = fullfile(folder_data, 'temp_wh.dat');

if ~isfile(path_data)
    error(['Required waveform file temp_wh.dat was not found in %s. ', ...
        'This pipeline currently requires temp_wh.dat and does not ', ...
        'automatically read raw .ap.bin files.'], folder_data);
end

file_info = dir(path_data);
if isempty(file_info) || file_info.bytes == 0
    error('Required waveform file temp_wh.dat is empty or unreadable in %s.', folder_data);
end

end
