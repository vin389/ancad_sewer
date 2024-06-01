function [t_idx, ht_m, rain_mm] = readHeightAndRain(filename)
    fprintf('# Reading data from file %s ...\n', filename)
    table = readtable(filename, 'Delimiter', ',', 'HeaderLines', 1);
    timeTxt = table.Var1;
    flow_m3_s = table.Var2;
    ht_m = table.Var3;
    rain_mm = table.Var4;
    nStep = size(table, 1);
    t_idx = linspace(1, nStep, nStep);
    t_days = (t_idx - 1) / 1440.; 
    fprintf('# Duration of data is %d minutes (%.1f days) starting from %s\n', ...
        nStep, nStep / 1440., datestr(timeTxt(1)));
end        
