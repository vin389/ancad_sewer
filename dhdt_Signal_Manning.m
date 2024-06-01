function dhdt = dhdt_Signal_Manning(gamma, lambda, hb, diameter, knss, rain_mm, ht_m)
% This function calculates change of water level (height) (per minute).

% Inputs:
%   gamma: coefficient of input (of rain)
%   lambda: coefficient of output (according to ht_m and hb when full tube)
%   hb: 
%   diameter: diameter of tube
%   knss: k/n*sqrt(ss) for manning equation.
%         k: 1.0 in SI unit
%         n: roughness coefficient in Manning equation
%         ss: stream slope
%   rain_mm: rain history (mm per .. hour?) (time series by minute) 
%   ht_m: water level (height) (time series by minute)
%
% Output:
%   dhdt: dhdt is a vector which is the same length with rain_mm and ht_m
    
    nt = size(rain_mm(:), 1); % number of time steps
    nt_ht = size(ht_m(:), 1); % number of time steps of height
    if nt ~= nt_ht
      error('dhdt_v2: rain length is %d. ht length is %d. They should be the same.', nt, nt_ht);
    end 
    % calculate the term related to water-in
    w_in = gamma * rain_mm(:); % w_in is a time series
    % calculate the term related to water-out
    w_out = zeros(nt, 1); 
    for i = 1: nt
        if (ht_m(i) >= diameter)
            % full tube
            w_out(i) = lambda * sqrt(max(0, ht_m(i) - hb));
        elseif (ht_m(i) < 0.0)
            % ht < 0, considered to be dry
            w_out(i) = 0.0;
        else
            % manning equation 
            [rh, af, pw] = circularFlowHydraulicRadius(0.5 * diameter, ht_m(i));
            vf = knss * rh ^ (2./3.);
            w_out(i) = af * vf; 
        end
    end
    % dhdt is water-in minus water-out
    dhdt = w_in - w_out; 
end
