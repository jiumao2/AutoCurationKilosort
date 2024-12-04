function plotClusterComparison(folder_data, spike_id_A, spike_id_B, title_name, save_filename)

    n_channel_plot = 4;
    n_random_spikes = 100;
    waveform_window = [-31,32];

    % load the data
    path_data = fullfile(folder_data, 'temp_wh.dat');
    amplitudes = readNPY(fullfile(folder_data, 'amplitudes.npy'));
    spike_clusters = readNPY(fullfile(folder_data, 'spike_clusters.npy'));
    spike_times = readNPY(fullfile(folder_data, 'spike_times.npy'));
    spike_templates = readNPY(fullfile(folder_data, 'spike_templates.npy'));
    load(fullfile(folder_data, 'ops.mat'));
    dir_output = dir(path_data);
    nFileSamp = dir_output.bytes ./ 2 ./ ops.Nchan;
    mmap = memmapfile(path_data, 'Format', {'int16', [ops.Nchan, nFileSamp], 'x'});

    spike_id_all = {spike_id_A, spike_id_B};
    waveforms_all = cell(1,2);
    mean_waveforms = cell(1,2);
    sorted_channels = cell(1,2);
    for k = 1:2
        spike_ids = spike_id_all{k};
        spike_times_this = spike_times(spike_ids);
    
        n_waveforms = min(length(spike_times_this), n_random_spikes);
        idx_rand = randperm(length(spike_times_this), n_waveforms);
        waveforms = zeros(n_waveforms, ops.Nchan, diff(waveform_window)+1); % nSpikes x 383 x 64
        for j = 1:n_waveforms
            waveforms(j,:,:) = mmap.Data.x(:,...
                spike_times_this(idx_rand(j)) + waveform_window(1):spike_times_this(idx_rand(j)) + waveform_window(2));
        end
        
        waveforms_all{k} = waveforms;
        mean_waveforms{k} = squeeze(mean(waveforms, 1)); % 383 x 64
        ptp = max(mean_waveforms{k}, [], 2) - min(mean_waveforms{k}, [], 2);
        [~, sorted_channels{k}] = sort(ptp, 'descend');
        
    end
    
    ch_included = [];
    for k = 1:ops.Nchan
        for j = 1:2
            ch_this = sorted_channels{j}(k);
            if ~any(ch_included == ch_this)
                ch_included = [ch_included, ch_this];
            end

            if length(ch_included) == n_channel_plot
                break
            end
        end
        if length(ch_included) == n_channel_plot
            break
        end
    end
    
    % compute the included channel ids
    channels = [1:ops.Nchan, 1:ops.Nchan];
    ptp = [max(mean_waveforms{1}, [], 2) - min(mean_waveforms{1}, [], 2);...
        max(mean_waveforms{2}, [], 2) - min(mean_waveforms{2}, [], 2)];
    
    [~, idx_sort] = sort(ptp, 'descend');
    ch_included = [];
    for k = 1:length(idx_sort)
        ch_this = channels(idx_sort(k));
        if ~any(ch_included == ch_this)
            ch_included = [ch_included, ch_this];
        end

        if length(ch_included) == n_channel_plot
            break
        end
    end
    

    fig = EasyPlot.figure();
    
    ax_waveform = EasyPlot.createGridAxes(fig, 1, n_channel_plot,...
        'Width', 3,...
        'Height', 3,...
        'MarginLeft', 0.1,...
        'MarginRight', 0.1,...
        'MarginBottom', 0,...
        'XAxisVisible', 'off',...
        'YAxisVisible', 'off');
    ax_autocorrelogram = EasyPlot.createGridAxes(fig, 2, 2,...
        'Width', 2,...
        'Height', 2,...
        'MarginLeft', 0.1,...
        'MarginRight', 0.1,...
        'MarginTop', 0.1,...
        'MarginBottom', 0.1,...
        'YAxisVisible', 'off');

    ax_amplitude = EasyPlot.axes(fig,...
        'Width', 8,...
        'Height', 3,...
        'MarginBottom', 1,...
        'YAxisVisible', 'off');

    EasyPlot.align(ax_autocorrelogram, ax_waveform, 'left');
    EasyPlot.place(ax_autocorrelogram, ax_waveform, 'bottom');
    EasyPlot.align(ax_amplitude, ax_waveform, 'right');
    EasyPlot.align(ax_amplitude, ax_autocorrelogram, 'verticalCenter');
    

    % overlapping waveforms on several nearest channels
    for k = 1:length(ax_waveform)
        plot(ax_waveform{k},...
            waveform_window(1):waveform_window(2),...
            squeeze(waveforms_all{1}(:,ch_included(k),:)), 'b-');
        plot(ax_waveform{k},...
            waveform_window(1):waveform_window(2),...
            squeeze(waveforms_all{2}(:,ch_included(k),:)), 'r-');
        title(ax_waveform{k}, ['Ch ', num2str(ch_included(k))]);
        xlim(ax_waveform{k}, [waveform_window(1), waveform_window(2)]);
    end
    EasyPlot.setYLim(ax_waveform);

    % cross correlogram
    spike_times_all = {spike_times(spike_id_A), spike_times(spike_id_B)};
    s = cell(1,2);
    max_t = max(max(spike_times_all{1} ./ 30), max(spike_times_all{2} ./ 30));
    for k = 1:length(s)
        temp = zeros(1, max_t);
        temp(spike_times_all{k} ./ 30) = 1;
        s{k} = temp;
    end

    for k = 1:2
        for j = 1:2
            [r, lags] = xcorr(s{k}, s{j}, 50);

            if k == j
                r(lags == 0) = 0;
            end
            
            if k == 1 && j == 1
                color_this = 'b';
            elseif k == 2 && j == 2
                color_this = 'r';
            else
                color_this = 'k';
            end
            bar(ax_autocorrelogram{k,j}, lags, r, 1, 'FaceColor', color_this, 'EdgeColor', 'none');
        end
    end

    EasyPlot.setXLim(ax_autocorrelogram, [-50, 50]);
    EasyPlot.set(ax_autocorrelogram, 'XAxisVisible', 'off');
    EasyPlot.set(ax_autocorrelogram{end,1}, 'XAxisVisible', 'on');

    % amplitude vs time
    min_t = min(min(spike_times_all{1}), min(spike_times_all{2}));
    max_t = max(max(spike_times_all{1}), max(spike_times_all{2}));
    amp_all = {amplitudes(spike_id_A), amplitudes(spike_id_B)};
    plot(ax_amplitude, double(spike_times_all{1}) ./ 30000, amp_all{1}, 'b.');
    plot(ax_amplitude, double(spike_times_all{2}) ./ 30000, amp_all{2}, 'r.');
    xlim(ax_amplitude, double([min_t, max_t])./30000);
    xlabel(ax_amplitude, 'Time (sec)');
    title(ax_amplitude, 'Amplitude')

    % rasters and PETH (if given)
    if exist(fullfile(folder_data, 'events.csv'), 'file') && exist(fullfile(folder_data, 'event_labels.csv'), 'file')
        n_max_event = n_channel_plot;
        t_pre = -1000;
        t_post = 1000;


        spike_times_all_ms = {double(spike_times(spike_id_A))./30,...
            double(spike_times(spike_id_B))./30};
        event_labels = readcell(fullfile(folder_data, 'event_labels.csv'));
        event_data = readcell(fullfile(folder_data, 'events.csv'));
        event_times_sec = cell(1, length(event_labels));
        n_event = min(n_max_event, length(event_labels));

        ax_raster = EasyPlot.createGridAxes(fig, 1, n_event,...
            'YAxisVisible', 'off',...
            'XAxisVisible', 'off',...
            'Width', 3,...
            'Height', 3,...
            'MarginLeft', 0.1,...
            'MarginRight', 0.1,...
            'MarginBottom', 1,...
            'MarginTop', 1);
        EasyPlot.align(ax_raster, ax_waveform, 'left');
        EasyPlot.place(ax_raster, ax_autocorrelogram, 'bottom');

        for k = 1:n_event
            for j = 1:size(event_data, 2)
                data_this = event_data{k,j};
                if ~ismissing(data_this)
                    event_times_sec{k} = [event_times_sec{k}, data_this];
                end
            end
        end

        x_plot = cell(2,n_event);
        y_plot = cell(2,n_event);
        for k = 1:n_event
            t_event_ms = event_times_sec{k}*1000;
            for j = 1:2
                st = spike_times_all_ms{j};
                for i = 1:length(t_event_ms)
                    t_event_this = t_event_ms(i);
                    st_this = st(st>t_event_this+t_pre & st<t_event_this+t_post)-t_event_this;
                    for ii = 1:length(st_this)
                        x_plot{j,k} = [x_plot{j,k}, st_this(ii), st_this(ii), NaN];
                        y_plot{j,k} = [y_plot{j,k}, i-0.5, i+0.5, NaN];
                    end
                end
            end
        end
        
        colors_all = {'b', 'r'};
        for k = 1:n_event
            for j = 1:2
                plot(ax_raster{k}, x_plot{j,k}, y_plot{j,k}, '-',...
                    'LineWidth', 1,...
                    'Color', colors_all{j});
            end

            ylim(ax_raster{k}, [0.5, length(event_times_sec{k})+0.5]);
            title(ax_raster{k}, event_labels{k});
        end

        EasyPlot.setXLim(ax_raster, [t_pre, t_post]);
        EasyPlot.set(ax_raster{1}, 'XAxisVisible', 'on');
    end

    
    h = EasyPlot.setGeneralTitle(ax_waveform, title_name);
    EasyPlot.move(h, 'dy', 0.3);

    output_folder = fileparts(save_filename);
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    EasyPlot.cropFigure(fig);
    EasyPlot.exportFigure(fig, save_filename, 'dpi', 600);
    close all;
end
