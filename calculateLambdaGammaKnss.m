function [h_b, lambda, gamma, knss, ht_t78_predict, dhdt_t78, ...
    dhdt_t78_predict] = ...
   calculateLambdaGammaKnss(t7, t8, rain_mm, ht_m, diameter, lambda_init, ...
   h_b_init, gamma_init, knss_init)
    % dhdt_t78
    t78 = double((t7:t8) - t7); t78 = t78(:);
    ht_t78 = ht_m(t7:t8); ht_t78 = ht_t78(:);   
    rain_t78 = rain_mm(t7:t8); rain_t78 = rain_t78(:);
    % calculate derivative (dh/dt) (diffy is dh/dt of selected time range)
    dhdt_t78_ = diff(ht_t78); dhdt_t78_ = dhdt_t78_(:);
    t78_diff_ = linspace(.5, double(t8 - t7) -.5, size(dhdt_t78_(:),1)); t78_diff_ = t78_diff_(:);
    dhdt_t78 = interp1(t78_diff_, dhdt_t78_, t78, 'linear', 'extrap');
    % define cost function 
    costFun = @(x) (dhdt_t78(:) - dhdt_Signal_Manning( ...
        x(1), x(2), x(3), diameter, x(4), rain_t78, ht_t78));
    % define initial guess, upper/lower bounds
    x0 = [gamma_init, lambda_init, h_b_init, knss_init];
    lb = [0.0, 0.0, -10.0, 0.0];
    ub = [1e6, 1e6,  10.0, 1e3];
    % run optimization (least squares)
    [x, resnorm, residual] = lsqnonlin(costFun, x0, lb, ub);
    % optimized parameters
    h_b = x(3);
    lambda = x(2);
    gamma = x(1);
    knss = x(4);
    % predict ht (t7 to t8)
    ht_t78_predict = ...
        predictHeight_Signal_Manning_ode23(t7, t8, rain_mm, ht_m, ...
        lambda, h_b, gamma, diameter, knss);
    dhdt_t78_predict = dhdt_Signal_Manning(...
        gamma, lambda, h_b, diameter, knss, ...
        rain_mm(t7:t8), ht_t78_predict);
end