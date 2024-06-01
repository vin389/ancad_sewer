function [t_idx, ht_m, rain_mm] = readHeightAndRain_smartiout_v7(filename)
    fprintf('# Reading data from single file (4-column, time/height/station/rain)\n')
    fprintf('#   from %s\n', filename);
    % read time data 
    table = readtable(filename, 'Delimiter', ',', 'HeaderLines', 1);
    timeTxt = table.Var1;
    nStep = size(timeTxt, 1);
    % read height data
    ht_m_txt = table.Var2;
    nStep_ht_m = size(ht_m_txt, 1);
    if (nStep_ht_m ~= nStep)
        fprintf('# Warning: Columns 1 and 2 have different # of rows (%d and %d). You could see some errors soon.\n', nStep, nStep_ht_m);
    end
    ht_m = zeros(nStep, 1);
    for i = 1: nStep
        strtmp = ht_m_txt{i};
        strtmp = strrep(strtmp,'(','');
        strtmp = strrep(strtmp,')','');
        strtmp = strrep(strtmp,'m','');
        ht_m(i) = str2double(strtmp);
    end
    % read rain data
    rain_mm_txt = table.Var4;
    nStep_rain_mm = size(rain_mm_txt, 1);
    if (nStep_rain_mm ~= nStep)
        fprintf('# Warning: Columns 1 and 4 have different # of rows (%d and %d). You could see some errors soon.\n', nStep, nStep_rain_mm);
    end
    rain_mm = zeros(nStep, 1);
    for i = 1: nStep
        strtmp = rain_mm_txt{i};
        strtmp = strrep(strtmp,' ','');
        strtmp = strrep(strtmp,'\t','');
        strtmp = strrep(strtmp,'(','');
        strtmp = strrep(strtmp,')','');
        strtmp = strrep(strtmp,'m','');
        rain_mm(i) = str2double(strtmp);
        if isnan(rain_mm(i)) 
            rain_mm(i) = 0.0;
        end
    end
    % for every 10 records, calculates average of rain
    for i = 11:10:nStep
        avgRain = sum(rain_mm(i-9:i)) / 10;
        rain_mm(i-9:i) = avgRain;
    end    
    % 
    t_idx = linspace(1, nStep, nStep);
end
