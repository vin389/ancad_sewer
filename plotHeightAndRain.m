function results = plotHeightAndRain(t_idx, ht_m, rain_mm, theAxes)
    % clear axes data 
    cla(theAxes); 
    % plot height (water level)
    axes(theAxes);
    plot(theAxes, t_idx, ht_m, 'b-', 'LineWidth', 2);  % Plot data2 on left y-axis
    xlabel(theAxes, 'Minutes');
    yyaxis(theAxes, 'left'); 
    ylabel(theAxes, 'Water level (height) (m)', 'Color', 'b');  % Label left y-axis in red
    % plot rainfall in the same axes. vertical label at right side
    hold(theAxes, 'on');
    plot(theAxes, t_idx, rain_mm, 'r--', 'LineWidth', 2);  % Plot rainfall on right y-axis
    yyaxis(theAxes, 'right');  % Activate right y-axis for subsequent plots
    ylabel(theAxes, 'Rainfall (mm)', 'Color', 'r');  % Label right y-axis
    title(theAxes, 'Rainfall and water level (height) vs time (day)');
    legend(theAxes, 'Water level (height)', 'Rainfall');
    grid(theAxes, 'on');
    xlim(theAxes, 'auto');
end
