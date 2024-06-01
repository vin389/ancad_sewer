function [gamma, cqin, lambda, hb, ckns, ...
    ht_sub_predict, dhdt_sub_diff, dhdt_sub_predict] = ...
   calculate_v3_qin_qout_params(rain_mm_sub, ht_m_sub, diameter, ...
    gamma_init, cqin_init, lambda_init, hb_init, ckns_init)
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
    costFun = @(x) (dhdt_sub_diff - dhdt_v3_qin_qout(...
        ht_m_sub, rain_mm_sub, diameter, ...
        x(1), x(2), x(3), x(4), x(5)) );
    % define initial guess, upper/lower bounds
    x0 = [gamma_init, cqin_init, lambda_init, hb_init, ckns_init];
%    lb = [0.0, 0.0,   0.0, -20.0, 0.0];
    hb_lower_bound = -100; 
    lb = [0.0, 0.0,   0.0, hb_lower_bound, 0.0];
    ub = [1e6, 1e6,  10.0,   0.0, 1e3];
    % run optimization (least squares)
    [x, resnorm, residual] = lsqnonlin(costFun, x0, lb, ub);
    % optimized parameters
    gamma = x(1);
    cqin = x(2);
    lambda = x(3);
    hb = x(4);
    ckns = x(5);
    % based on the optimized parameters
    %   predict ht (ht_m_sub, the ht_m between selected time range)
    ht_init = ht_m_sub(1);
    ht_sub_predict = ...
        predictHeight_v3_qin_qout(...
            rain_mm_sub, ht_init, ...
            diameter, gamma, cqin, lambda, hb, ckns); 
    % based on the optimized parameters
    %   based on the physical model, calculate the dhdt 
    dhdt_sub_predict = dhdt_v3_qin_qout(...
        ht_m_sub, rain_mm_sub, diameter, ...
        gamma, cqin, lambda, hb, ckns);
end