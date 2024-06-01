function dhdt = dhdt_v4_p2_qout(h0, diameter, alpha)
% This function calculates change of water level (height) (per minute).
% Model:
%   Q_inlet = 0
%   Q_outlet
%     h > D: Q_outlet = -alpha * h0 
%     h <= D: 0
%   if h > D: (assuming (dh/dt) is proportional to Q)
%     dhdt = C_0 * (Q_inlet - Q_outlet) 
%          = -(C_0 * -alpha' * h0)
%          = alpha * h0          
%   if h <= D: (assuming Q is proportional to h^1.68) (?)
%              or (assuming (dh/dt) is proportional to Q )
%     dhdt = 0
%     

% Inputs:
%   h0: main part of ht water level (height) (time series by minute)
%   diameter: diameter of tube
%   alpha: alpha = -0.5 * lambda / sqrt(-hb)
% Output:
%   dhdt: dhdt is a vector which is the same length with h0
%     Note: dhdt is based on the model. It is not calculated by numrical
%           derivative like finite difference. 
%           This function is used to tune (by optimization) parameters like
%           alpha
%    
    nt = size(h0(:), 1); % number of time steps of height
    
    % logic arrays. isFullTube(i) is true (1) if ht_m(i) > diameter.
    %               notFullTube(i) is inverse of isFullTube
    %               This prevents an explicit for-loop in our later matlab code. 
    isFullTube = h0(:) >= diameter;
    % calculate the term related to water-out
    %   v3: w_out_full = lambda * sqrt(max(0, ht_m(:) - hb));
    %   v4 p2: w_out_full = -alpha * h0;
    w_out_full = -alpha * h0;
    % combine two cases
    dhdt = isFullTube  .* (- w_out_full);
end
