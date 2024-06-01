% This function is only for function predictHeight_ode23. 
% Function predictHeight_ode23 calls Matlab ode23, which requires a 
% defined function of odefun.
% ode23 requires an ode function (odefun) that is defined by user
% Note: odefun is a Matlab function, not an array, but in our case, 
%       the function contains rain-fall, which is an array, not a function,
%       so we need to use interp1 to convert array rain_t56 to a value of 
%       a function, especially while the t could be not an integer.
function dhdt = predictHeight_ode23_odefun(...
                t, ht_t, t56, rain_t56, lambda, hb, gamma)
    % while t can be a float, we need to interpolate values in rain_t56
    rain_t = interp1(t56, rain_t56, t); 
    % calculate ode function (, i.e., y' = ode function)
    dhdt = gamma * rain_t - lambda * ... 
           sqrt(ht_t - hb);
end
