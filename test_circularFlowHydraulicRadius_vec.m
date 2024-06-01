function test_circularFlowHydraulicRadius_vec()
    clear rks; clear afs; clear pws; clear depths; 
    % define data
    tic;
    r = 100.;
    N = 10000000; % number of points to plot
    depths=linspace(0, 2 * r, N); 
    depths = depths(:);
    [rks, afs, pws] = circularFlowHydraulicRadius_vec(r, depths);
    tcost = toc;
    fprintf('# Computing time of vector version: %f sec.\n', tcost);
    % Create a figure and arrange subplots
    figure;

    % Subplot 1
    subplot(3, 1, 1);  % 3 rows, 1 column, subplot 1
    plot(depths/r, rks / r);
    title('Hydraulic radius');
    xlabel('depth / r');
    ylabel('Rk / r'); 
    grid('on');

    % Subplot 2
    subplot(3, 1, 2);  % 3 rows, 1 column, subplot 2
    plot(depths/r, afs/(pi*r*r));
    title('Flow area');
    xlabel('depth / r');
    ylabel('Af / (pi*r*r)');
    grid('on');

    % Subplot 3 
    subplot(3, 1, 3);  % 3 rows, 1 column, subplot 3
    plot(depths / r, pws / r);
    title('Wetted perimeter');
    xlabel('depth/r');
    ylabel('pw/r');
    grid('on');

    % Adjust layout (optional)
    sgtitle('Checking hydraulic Radius (Vec version)');  % Add a mai
end