
function dhdt = dhdt_v3_qin_qout_ode23(...
                t, ht_t, rain_mm_sub, ...
                diameter, gamma, cqin, lambda, hb, ckns)
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
%   diameter, gamma, cqin, lambda, hb, ckns: parameters: optimized 
%      parameters.  
    % while t can be a float, we need to interpolate values in rain_mm_sub
    nt = size(rain_mm_sub(:), 1);
    rain_t = interp1(t910, rain_t910, t); 
    % calculate ode function (, i.e., y' = ode function)
    dhdt = dhdt_Signal_Manning(gamma, lambda, hb, diameter, knss, rain_t, ht_t);
end
