function dhdt = dhdt_v4_p1_qin_qout(ht_m, rain_mm, diameter, gamma, cqin, kappa, ckns)
% This function calculates change of water level (height) (per minute).
% Model:
%   Q_inlet = gamma' * rain(t) + Q_inNormal (or Q_inNormal(t))
%   Q_outlet
%     v4: h > D: Q_outlet = kappa 
%       (compare to v3: h > D: Q_outlet = lambda' * sqrt(h - hb)  )
%     h <= D: Q_outlet = k/n * af * rh^(2/3) * sqrt(ss)
%                      = knss * af * rh^(2/3)
%   if h > D: (assuming (dh/dt) is proportional to Q)
%     dhdt =  C_0 * (Q_inlet - Q_outlet) 
%          = + (C_0 * gamma') rain(t) + C_0 * Q_inNormal
%            - (C_0 * kappa')
%          =  gamma * rain(t) + cqin 
%           - kappa           
%   if h <= D: (assuming Q is proportional to h^1.68) (?)
%              or (assuming (dh/dt) is proportional to Q )
%     dhdt =  gamma * rain(t) + cqin 
%           - ckns * af * rh^(2/3)
%       where
%          cqin is a parameter, which is from C_0 * Q_inNormal
%          ckns is a parameter, which is from C_0 * k/n * sqrt(stream slope)
%     
%     Flow area (af) and hydrolic radius (rh) can be calculated by 
%     the function:
%      [rh, af, pw] = circularFlowHydraulicRadius(r, depth)
%        where r is radius of (horizontal) tube
%              depth is about the same as h(t), the water level 
%  Manning equation: 
%    Q = VA = 1/n * A * R^(2/3) * sqrt(stream slope)
%
 
% Inputs:
%   ht_m: water level (height) (time series by minute)
%   rain_mm: rain history (mm per .. hour?) (time series by minute) 
%   diameter: diameter of tube
%   gamma: coefficient of input (of rain)
%   cqin: a parameter related to Q_inlet (regular sewage) (C_0 * Q_inlet)
%   kappa: coefficient of output. kappa is close to lambda*sqrt(-hb)
%   ckns: a simplified parameter related to Manning equation, which is
%         C_0 * k/n * sqrt(ss), where
%         k: 1.0 in SI unit
%         n: roughness coefficient in Manning equation
%         ss: stream slope
%
% Output:
%   dhdt: dhdt is a vector which is the same length with rain_mm and ht_m
%     Note: dhdt is based on the model. It is not calculated by numrical
%           derivative like finite difference. 
%           This function is used to tune (by optimization) parameters like
%           gamma, cqin, lambda, hb, and ckns
%    
    nt = size(rain_mm(:), 1); % number of time steps
    nt_ht = size(ht_m(:), 1); % number of time steps of height
    if nt ~= nt_ht
      error('dhdt_v4_p1_qin_qout: rain length is %d. ht length is %d. They should be the same.', nt, nt_ht);
    end 
    
    % logic arrays. isFullTube(i) is true (1) if ht_m(i) > diameter.
    %               notFullTube(i) is inverse of isFullTube
    %               This prevents an explicit for-loop in our later matlab code. 
    isFullTube = ht_m(:) >= diameter;
    notFullTube = ~isFullTube;
    % calculate the term related to water-in
    w_in_full = gamma * rain_mm(:) + cqin; % w_in is a time series
    w_in_notf = w_in_full; 
    % calculate the term related to water-out
    %   v3: w_out_full = lambda * sqrt(max(0, ht_m(:) - hb));
    %   v4: w_out_full = kappa;
    w_out_full = kappa;
    [rhs, afs, pws] = circularFlowHydraulicRadius_vec(diameter / 2., ht_m(:));
    vfs = ckns * rhs(:) .^ (2./3.);
    w_out_notf = afs(:) .* vfs(:); 
    % combine two cases
    dhdt = isFullTube  .* (w_in_full - w_out_full) ...
         + notFullTube .* (w_in_notf - w_out_notf);
end
