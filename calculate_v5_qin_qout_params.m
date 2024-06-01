function [hb, kappa, gamma, ckns, cqin, ...
    ht_sub_predict, dhdt_sub_diff, dhdt_sub_predict] = ...
   calculate_v5_qin_qout_params(ht_m_sub, rain_mm_sub, diameter, ...
    hb_init, kappa_init, gamma_init, ckns_init, cqin_init)
%    % dhdt_t78
%    ht_t78 = ht_m(t7:t8); ht_t78 = ht_t78(:);   
%    rain_t78 = rain_mm(t7:t8); rain_t78 = rain_t78(:);
    % calculate derivative (dh/dt) (diffy is dh/dt of selected time range)
    nt = size(ht_m_sub(:), 1); 
    dhdt_sub_diff_ = diff(ht_m_sub); dhdt_sub_diff_ = dhdt_sub_diff_(:);
    dhdt_sub_t_ = linspace(1.5, nt - 0.5, nt - 1); dhdt_sub_t_ = dhdt_sub_t_(:);
    dhdt_sub_diff = interp1(dhdt_sub_t_, dhdt_sub_diff_, 1:nt, 'linear', 'extrap');
    dhdt_sub_diff = dhdt_sub_diff(:);
    % define cost function, which calls dhdt_xxxx (the physical model)
    costFun = @(x) (dhdt_sub_diff - dhdt_v5_qin_qout(...
        ht_m_sub, rain_mm_sub, diameter, ...
        x(1), x(2), x(3), x(4), x(5)) );
    % define initial guess, upper/lower bounds
    x0 = [hb_init, kappa_init, gamma_init, ckns_init, cqin_init];
    hb_lower_bound = -100.; 
    lb = [hb_lower_bound, 0, 0, 0 ,0];
    ub = [0.0, 1e6, 1e6, 1e6, 1e6, 1e6];
%    kappa_lb = 0.025; 
%    kappa_ub = 0.035; 
%    gamma_lb = 0.011;
%    gamma_ub = 0.011; 
%    lb = [hb_lower_bound, kappa_lb, gamma_lb, 0 ,0];
%    ub = [0.0,            kappa_ub, gamma_ub, 1e6, 1e6];
    % run optimization (least squares)
    [x, resnorm, residual] = lsqnonlin(costFun, x0, lb, ub);
    % optimized parameters
    hb = x(1);
    kappa = x(2);
    gamma = x(3);
    ckns = x(4);
    cqin = x(5);
    fprintf('The LSQ residual is %f\n', norm(residual));
    % based on the optimized parameters
    %   predict ht (ht_m_sub, the ht_m between selected time range)
    ht_init = ht_m_sub(1);
    ht_sub_predict = ...
        predictHeight_v5_qin_qout(...
            ht_init, rain_mm_sub, diameter, ...
            hb, kappa, gamma, ckns, cqin); 
    % based on the optimized parameters
    %   based on the physical model, calculate the dhdt 
    dhdt_sub_predict = dhdt_v5_qin_qout(...
        ht_m_sub, rain_mm_sub, diameter, ...
        hb, kappa, gamma, ckns, cqin);
end