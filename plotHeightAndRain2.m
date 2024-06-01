function results = plotHeightAndRain2(t_idx, ht_m, rain_mm, theAxes1, theAxes2)
    % clear axes data 
    cla(theAxes1); xlim(theAxes1, [0 1]); 
    cla(theAxes2); xlim(theAxes2, [0 1]); 
    % plot height (water level)
    plot(theAxes1, t_idx, ht_m, 'b', 'LineWidth', 2);  
    xlabel(theAxes1, 'Minutes');
    ylabel(theAxes1, 'Height (water level) (m)', 'Color', 'black');
    grid(theAxes1, 'on');
    xlim(theAxes1, 'auto');
    % plot rainfall 
    plot(theAxes2, t_idx, rain_mm, 'r', 'LineWidth', 2); 
    ylabel(theAxes2, 'Rainfall (mm)', 'Color', 'black');  
    grid(theAxes2, 'on');
    xlim(theAxes2, 'auto');
    % link axes of these plots
    linkaxes([theAxes1 theAxes2], 'x');
end