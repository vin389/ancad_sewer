function [t_idx, ht_m, rain_mm] = readHeightAndRain_twoFiles(source_dir, source_file)
    fprintf('# Reading data from two files two files (rain and sensors)\nDir:%s\nRain:%s\nSewer:%s ...\n', source_dir, source_file{1}, source_file{2})
    % read rain data 
    file1 = [source_dir source_file{1}];
    table1 = readtable(file1, 'Delimiter', ',', 'HeaderLines', 1);
    timeTxt1 = table1.Var1;
    time1 = zeros(length(timeTxt1),1);
    for i = 1:length(time1)
        % datenum uses unit of day. We need unit of minute. 
        time1(i) = datenum(timeTxt1(i)) * 1440; 
    end
    rain_mm_nonsync = table1.Var2;
    % read sensor data (height)
    file2 = [source_dir source_file{2}];
    table2 = readtable(file2, 'Delimiter', ',', 'HeaderLines', 1);
    timeTxt2 = table2.Var1;
    time2 = zeros(length(timeTxt2),1);
    for i = 1:length(time2)
        % datenum uses unit of day. We need unit of minute. 
        time2(i) = datenum(timeTxt2(i)) * 1440; 
    end
    ht_m_nonsync = table2.Var2;
    % synchronization 
    t_start = round(min([time1(1), time2(1)]));
    t_end = round(max([time1(end), time2(end)]));
    nStep = t_end - t_start + 1;
    t_sync = linspace(t_start, t_end, nStep);
    % check rain_mm_nonsync has nan. 
    % if yes, do correction by simple interpolation.
    for k = 1:size(rain_mm_nonsync(:), 1)
        if isnan(rain_mm_nonsync(k))
            fprintf('# Warning: rain data %d has an nan problem. Check your file.\n', k);
            if k > 1 && k < size(rain_mm_nonsync(:), 1)      
                time1(k) = .5 * (time1(k-1) + time1(k+1))
                rain_mm_nonsync(k) = .5 * (rain_mm_nonsync(k-1) + rain_mm_nonsync(k+1))
            elseif k == 1
                time1(k) = 2 * time1(k+1) - time1(k+2)
                rain_mm_nonsync(k) = 2 * (rain_mm_nonsync(k+1) - rain_mm_nonsync(k+2))
            elseif k == size(rain_mm_nonsync(:), 1)
                time1(k) = 2 * time1(k-1) - time1(k-2)
                rain_mm_nonsync(k) = 2 * (rain_mm_nonsync(k-1) - rain_mm_nonsync(k-2))
            end
        end
    end
    % matlab interpolation 
    rain_mm = interp1(time1, rain_mm_nonsync, t_sync, 'nearest', 'extrap');
    % check rain_mm_nonsync has nan. 
    % if yes, do correction by simple interpolation.
    for k = 1:size(ht_m_nonsync(:), 1)
        if isnan(ht_m_nonsync(k))
            fprintf('# Warning: water height data %d has an nan problem. Check your file.\n', k);
            if k > 1 && k < size(ht_m_nonsync(:), 1)      
                time2(k) = .5 * (time2(k-1) + time2(k+1))
                ht_m_nonsync(k) = .5 * (ht_m_nonsync(k-1) + ht_m_nonsync(k+1))
            elseif k == 1
                time2(k) = 2 * time2(k+1) - time2(k+2)
                ht_m_nonsync(k) = 2 * (ht_m_nonsync(k+1) - ht_m_nonsync(k+2))
            elseif k == size(rain_mm_nonsync(:), 1)
                time2(k) = 2 * time2(k-1) - time2(k-2)
                ht_m_nonsync(k) = 2 * (ht_m_nonsync(k-1) - ht_m_nonsync(k-2))
            end
        end
    end
    % matlab interpolation 
    ht_m = interp1(time2, ht_m_nonsync, t_sync, 'nearest', 'extrap');
    t_idx = linspace(1, nStep, nStep);
    % print info
    if (time1(1) < time2(1))
        t_start_txt = timeTxt1(1);
    else
        t_start_txt = timeTxt2(1);
    end
    fprintf('# Duration of data is %d minutes (%.1f days) starting from %s\n', ...
        nStep, nStep / 1440., datestr(t_start_txt));
end        
