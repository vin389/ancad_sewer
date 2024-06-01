function [gamma, cqin, lambda, hb, ckns, ...
    ht_sub_predict, dhdt_sub_diff, dhdt_sub_predict] = ...
   calculate_v4_qin_qout_params(rain_mm_sub, ht_m_sub, diameter, ...
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
    costFun = @(x) (dhdt_sub_diff - dhdt_v4_p1_qin_qout(...
        ht_m_sub, rain_mm_sub, diameter, ...
        x(1), x(2), x(3), x(4)) );
    % define initial guess, upper/lower bounds 
    %   note: if hb -> 0, alpha -> infinite, so we replace 
    %         sqrt(-hb_init) to sqrt(max(1, -hb_init))  
    kappa_init = lambda_init * sqrt(max(1., -hb_init));
    x0 = [gamma_init, cqin_init, kappa_init, ckns_init];
%    lb = [0.0, 0.0,   0.0, -20.0, 0.0];
    lb = [0.0, 0.0,  0.0,  0.0];
    ub = [1e6, 1e6,  1e6,  1e3];
    % run optimization (least squares)
    [x, resnorm, residual] = lsqnonlin(costFun, x0, lb, ub);
    % optimized parameters
    gamma = x(1);
    cqin = x(2);
    kappa = x(3);
    ckns = x(4);
    % based on the optimized parameters
    %   predict ht (ht_m_sub, the ht_m between selected time range)
    % in v4, the predicted ht is actually ht0
    ht_init = ht_m_sub(1);
    h0_sub_predict = ...
        predictHeight_v4_p1_qin_qout(...
            rain_mm_sub, ht_init, ...
            diameter, gamma, cqin, kappa, ckns); 
    %
    % Phase 2 (p2)
    % eh1 = h - h0
    %  (the _sub indicates the array is not entire history, instead, 
    %   it is only duration of t7 to t8)
    eh1_sub = ht_m_sub - h0_sub_predict; 
    % calculate derivative (deh1/dt) 
    nt = size(eh1_sub(:), 1); 
    deh1dt_sub_diff_ = diff(eh1_sub); deh1dt_sub_diff_ = deh1dt_sub_diff_(:);
    deh1dt_sub_t_ = linspace(1.5, nt - 0.5, nt - 1); deh1dt_sub_t_ = deh1dt_sub_t_(:);
    deh1dt_sub_diff = interp1(deh1dt_sub_t_, deh1dt_sub_diff_, 1:nt, 'linear', 'extrap');
    deh1dt_sub_diff = deh1dt_sub_diff(:);
    % define cost function, which calls dhdt_xxxx (the physical model)
    costFun_p2 = @(x) (deh1dt_sub_diff - dhdt_v4_p2_qout(...
        h0_sub_predict, diameter, x) );
    % define initial guess, upper/lower bounds 
    %   note: if hb -> 0, alpha -> infinite, so we replace 
    %         sqrt(-hb_init) to sqrt(max(1, -hb_init))  
    alpha_init = -0.5 * lambda_init / sqrt(max(1., -hb_init));
    x0 = alpha_init;
%    lb = [0.0, 0.0,   0.0, -20.0, 0.0];
    lb = -1e6;
    ub = -1e-6;
    % run optimization (least squares)
    [x, resnorm, residual] = lsqnonlin(costFun_p2, x0, lb, ub);
    alpha = x;
    % calculate lambda and hb
    lambda = sqrt(max(0, -2 * alpha * kappa));
    hb = -1. * (kappa/lambda)^2; 
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