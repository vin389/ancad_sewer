function [lambda, hb, dhdt_t12, lambda_corrcoef, ht_simu, t_hzero, p] = calculateLambda(t1, t2, ht_m)
   % copy data from water level
   t12 = double((t1:t2) - t1); t12 = t12(:);
   ht_t12 = ht_m(t1:t2); ht_t12 = ht_t12(:);
   % calculate derivative (dh/dt) (diffy is dh/dt of selected time range)
   dhdt_t12_ = diff(ht_t12); dhdt_t12_ = dhdt_t12_(:);
   t12_diff_ = linspace(.5, double(t2 - t1) -.5, size(dhdt_t12_(:),1)); t12_diff_ = t12_diff_(:);
   dhdt_t12 = interp1(t12_diff_, dhdt_t12_, t12, 'linear', 'extrap');
   % linear regression to get the slope of dh_dt, which is 0.5*(lambda^2)
   % and solve lambda and hb 
   p = polyfit(t12, dhdt_t12, 1);
   corrmat = corrcoef(t12, dhdt_t12);
   lambda_corrcoef = corrmat(1,2);
   lambda = sqrt(2 * p(1));
   h0 = ht_m(t1);
   hb = h0 - (-p(2)/lambda)^2;
   % simulated h(t), based on the linear regression result p()
   ht_simu = hb + (-0.5 * lambda * t12 + sqrt(h0 - hb)).^2;
   % time t_hzero is the time where h(t) is zero
   t_hzero = 2. / lambda * (sqrt(h0 - hb) - sqrt(-hb)); 
end