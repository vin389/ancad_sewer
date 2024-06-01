function [gamma, ht_t34_predict, dhdt_t34, dhdt_t34_fitting] = ...
   calculateGamma(t3, t4, rain_mm, ht_m, lambda, hb, h_min_1, h_min_2)
   t34 = double((t3:t4) - t3); t34 = t34(:);
   ht_t34 = ht_m(t3:t4); ht_t34 = ht_t34(:);   
   rain_t34 = rain_mm(t3:t4); rain_t34 = rain_t34(:);
   % calculate derivative (dh/dt) (diffy is dh/dt of selected time range)
   dhdt_t34_ = diff(ht_t34); dhdt_t34_ = dhdt_t34_(:);
   t34_diff_ = linspace(.5, double(t4 - t3) -.5, size(dhdt_t34_(:),1)); t34_diff_ = t34_diff_(:);
   dhdt_t34 = interp1(t34_diff_, dhdt_t34_, t34, 'linear', 'extrap');
   % define cost function for solving gamma 
   % That is, finding best gamma so that, 
   % dhdt = gamma * rain - lambda * sqrt(h - hb) 
   eq28_ver = 3;
   if eq28_ver == 1
      lambda_lost_term =  lambda * sqrt(max(0., ht_t34 - hb)); 
   elseif eq28_ver == 2
      lambda_lost_term =  lambda * sqrt(max(0., ht_t34 - hb)) .* heaviside(ht_t34 - 0.3); 
   elseif eq28_ver == 3
      lambda_lost_term =  lambda * sqrt(max(0., ht_t34 - hb)) .* piecewise_linear(h_min_1, h_min_2, ht_t34); 
   end
   err_dhdt = @(gamma) ... 
		   (dhdt_t34 - (gamma * rain_t34 - lambda_lost_term));
 
   % run nonlinear least squares
   gamma_init = 0.0;
   gamma_lb = 0.0; 
   gamma_ub = 1e6;
   [gamma,resnorm,residual,exitflag,output,lsq_lambda,jacobian] = ...
       lsqnonlin(err_dhdt, gamma_init, gamma_lb, gamma_ub);
   
   % debug for checking optimized gamma 
   theDebug = 0; 
   if (theDebug >= 1) 
       gamma_low = gamma * 0.95;
       io_low = (gamma_low * rain_t34 - lambda_lost_term);
       gamma_opt = gamma * 1.00;
       io_opt = (gamma_opt * rain_t34 - lambda_lost_term);
       gamma_hi = gamma * 1.05;
       io_hi = (gamma_hi * rain_t34 - lambda_lost_term);
       figure; plot(t34, dhdt_t34, t34, io_low, t34, io_opt, t34, io_hi); legend('dhdt', 'I/O low Gamma', 'I/O opt Gamma', 'I/O hi Gamma'); 
   end 
   
   % for debug: calculate residual 	   
%   residual_ = dhdt_t34 - (gamma * rain_t34 - lambda_lost_term);
%   norm_residual = norm(dhdt_t34 - (gamma * rain_t34 - lambda_lost_term));
%   figure; plot(t34 + double(t3), dhdt_t34, t34 + double(t3), (gamma * rain_t34 - lambda_lost_term)); legend('Measured', 'Fitting');
   % calculate predicted htdt_t34 (htdt_t34_predict) based on lambda, hb, gamma
   dhdt_t34_fitting = gamma * rain_t34 - lambda_lost_term;
   % calculate predicted ht_t34 (ht_t34_predict) based on lambda, hb, gamma
   ht_t34_predict = predictHeight_ode23(t3, t4, rain_mm, ht_m, lambda, hb, gamma)
   
end


	