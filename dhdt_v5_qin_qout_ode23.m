
function dhdt = dhdt_v5_qin_qout_ode23(...
                t, ht_t, rain_mm_sub, diameter, ...
                hb, kappa, gamma, ckns, cqin)
% This function is only for function predictHeight_Signal_Manning_ode23. 
% Function predictHeight_Signal_Manning_ode23 calls Matlab ode23, which 
% requires a defined function of odefun.
% ode23 requires an ode function (odefun) that is defined by user
% Note: odefun is a Matlab function, not an array, but in our case, 
%       the function contains rain-fall, which is an array, not a function,
%       so we need to use interp1 to convert array rain_t56 to a value of 
%       a function, especially while the t could be not an integer.            
%
% Input: 
%   t: the time. It is the same as index of array, only t is a float not an
%      integer. That means this function needs to do interpolation.
%      The reason that t is a float is, this function is called by ode23
%      and ode23 can input arbitrary value of t, not intger. 
%   ht_t: the water level at given time t, determined by ode23, as this
%      function is called by ode23. 
%   rain_mm_sub: the rain data. rain_mm_sub(1) is the rain data where t=1,
%      and so on.
%   diameter:
%   hb, kappa, gamma, ckns, cqin: parameters: optimized parameters.  
    % while t can be a float, we need to interpolate values in rain_mm_sub
    nt = size(rain_mm_sub(:), 1);
    rain_t = interp1(linspace(1, nt, nt), rain_mm_sub, t); 
    % calculate ode function (, i.e., y' = ode function)
    dhdt = dhdt_v5_qin_qout(ht_t, rain_t, diameter, ...
        hb, kappa, gamma, ckns, cqin);
end
